"""
s3_atlas.py -- read the atlas_of_chromatin from S3 for the NIPT (seqFF++) fork.

Design: METADATA-FIRST. Read the per-modality index files (encode_<mod>_info.txt)
once to learn what exists, filter in that table, then GET only the peak files that
survive the filter. No file in the atlas is opened until it is actually selected.

Auth: uses the default boto3 credential chain -- on EC2/SageMaker that is the
attached IAM role (no keys in code). Nothing here downloads patient data; it only
reads the public-within-org chromatin atlas.

Bucket / layout (confirmed):
  s3://<BUCKET>/niptai/data/atlas_of_chromatin/
      <modality>/ENCODE/bed_narrowPeak_hg19/ENCFF*.csv.gz      # peak files
      metadata/ENCODE/encode_<modality>_info.tsv               # the index (TSV,
      metadata/ENCODE/encode_<modality>_links.txt               #  1st line = report URL)
      ref/hg38ToHg19.over.chain
  modality in {chip, atac, DNase, scatac, SEdb, nucmap, dhs_index, sedb}

Peak-file columns (self-describing):
  seqnames,start,end,width,strand,name,score,signalValue,pValue,qValue,peak,
  Accession,Assay title,Target of assay,Biosample_class,Biosample_name,Genome
  -> use signalValue (float fold-enrichment), NOT score (saturated at 1000).
"""
from __future__ import annotations

import gzip
import io
import os
import re
import time

import pandas as pd

try:  # boto3 only needed for actual S3 access, not for importing the module
    from botocore.exceptions import ClientError
except ImportError:  # allows offline import of downstream pure-numeric code
    class ClientError(Exception):
        pass

BUCKET = os.environ.get("ATLAS_BUCKET", "dgx-bx-okapi-stage-production")
ATLAS_PREFIX = os.environ.get("ATLAS_PREFIX", "niptai/data/atlas_of_chromatin")

# columns we actually keep from each peak file
PEAK_USECOLS = ["seqnames", "start", "end", "signalValue",
                "Accession", "Assay title", "Target of assay",
                "Biosample_class", "Biosample_name", "Genome"]


def _client():
    # default chain: IAM role on EC2/SageMaker, or ~/.aws / env locally
    import boto3  # lazy: only required when S3 is actually touched
    return boto3.client("s3")


# ── metadata-first catalog ───────────────────────────────────────────────────
def read_index(modality: str = "chip", s3=None) -> pd.DataFrame:
    """Read metadata/ENCODE/encode_<modality>_info.tsv into a DataFrame.

    This is the ONE small read per modality that replaces scanning every peak file.
    Returns the experiment-level index (one row per ENCODE experiment, ENCSR*).

    The ENCODE dumps have TWO quirks handled here:
      * the first line is the ENCODE /report/ query URL, not the header -- the real
        tab-delimited header (ID, Accession, Assay name, ...) is on the SECOND line;
      * the file is tab-separated and named .tsv (older dumps used .txt).
    """
    s3 = s3 or _client()
    body = None
    for ext in (".tsv", ".txt"):
        key = f"{ATLAS_PREFIX}/metadata/ENCODE/encode_{modality}_info{ext}"
        try:
            body = s3.get_object(Bucket=BUCKET, Key=key)["Body"].read()
            break
        except ClientError as e:
            if e.response["Error"]["Code"] in ("NoSuchKey", "404"):
                continue
            raise
    if body is None:
        raise FileNotFoundError(f"no encode_{modality}_info.[tsv|txt] under {ATLAS_PREFIX}/metadata/ENCODE/")

    txt = body.decode("utf-8", errors="replace")
    lines = txt.split("\n")
    # find the header row: first line containing the 'Accession' column token
    hdr = next((i for i, ln in enumerate(lines[:5]) if "Accession" in ln and "\t" in ln), 0)
    return pd.read_csv(io.StringIO("\n".join(lines[hdr:])), sep="\t", dtype=str)


_ENCFF_RE = re.compile(r"ENCFF[0-9A-Z]+")


def list_bed_files(modality: str = "chip", genome: str = "hg19", s3=None) -> set:
    """List the ENCFF file accessions physically present in the hg19 bed folder.

    One paginated list_objects (metadata-only, no GET/decompress) -- this is how we
    learn which of an experiment's many files is the processed narrowPeak bed,
    instead of scanning every file's header.
    """
    s3 = s3 or _client()
    prefix = f"{ATLAS_PREFIX}/{modality}/ENCODE/bed_narrowPeak_{genome}/"
    accs = set()
    token = None
    while True:
        kw = dict(Bucket=BUCKET, Prefix=prefix)
        if token:
            kw["ContinuationToken"] = token
        resp = s3.list_objects_v2(**kw)
        for obj in resp.get("Contents", []):
            m = _ENCFF_RE.search(obj["Key"].rsplit("/", 1)[-1])
            if m:
                accs.add(m.group(0))
        if not resp.get("IsTruncated"):
            break
        token = resp.get("NextContinuationToken")
    return accs


def build_catalog(modalities=("chip", "atac", "DNase"), genome="hg19",
                  s3=None) -> pd.DataFrame:
    """Union the per-modality indexes into one FILE-level catalog for `genome`.

    The metadata index is EXPERIMENT-level (one ENCSR* row lists many ENCFF* files
    in the `Files` column, across assemblies/types). The peak files are named by
    FILE accession (ENCFF*). So we:
      1. read the experiment index (read_index),
      2. list the ENCFF accessions actually present in the hg19 bed folder
         (list_bed_files -- one metadata list, no header scanning),
      3. explode `Files` to ENCFF ids and keep only those present in the folder.

    Each surviving row carries the experiment-level annotation (tissue, target,
    assay) AND a resolved `s3_key` to the actual bed. Cache to parquet so this
    runs once, never per analysis.
    """
    s3 = s3 or _client()
    frames = []
    for mod in modalities:
        try:
            idx = read_index(mod, s3=s3)
        except (ClientError, FileNotFoundError) as e:
            print(f"[catalog] skip {mod}: {e}")
            continue
        try:
            present = list_bed_files(mod, genome=genome, s3=s3)
        except ClientError as e:
            print(f"[catalog] skip {mod} (no bed folder): {e.response['Error']['Code']}")
            continue

        idx = idx.copy()
        idx.columns = [c.strip() for c in idx.columns]
        # optional experiment-level assembly filter (rows list e.g. 'GRCh38,hg19')
        gcol = next((c for c in idx.columns if c.lower() in ("genome assembly", "genome", "assembly")), None)
        if gcol and genome:
            idx = idx[idx[gcol].astype(str).str.contains(genome, case=False, na=False)]

        files_col = next((c for c in idx.columns if c.strip().lower() == "files"), None)
        exp_col = next((c for c in idx.columns if c.strip().lower() == "accession"), None)
        if files_col is None or exp_col is None:
            print(f"[catalog] skip {mod}: no Files/Accession column ({list(idx.columns)[:6]}...)")
            continue

        rows = []
        annot_cols = [c for c in idx.columns if c != files_col]
        for _, r in idx.iterrows():
            for facc in set(_ENCFF_RE.findall(str(r[files_col]))):
                if facc in present:
                    rec = {c: r[c] for c in annot_cols}
                    rec["file_accession"] = facc
                    rec["experiment"] = r[exp_col]
                    rec["modality"] = mod
                    rec["s3_key"] = (f"{ATLAS_PREFIX}/{mod}/ENCODE/"
                                     f"bed_narrowPeak_{genome}/{facc}_bed.csv.gz")
                    rows.append(rec)
        frames.append(pd.DataFrame(rows))
        print(f"[catalog] {mod}: {len(idx)} experiments -> {len(rows)} hg19 bed files")

    cat = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    print(f"[catalog] {len(cat)} total bed files across {list(modalities)} ({genome})")
    return cat


# ── archived-tier handling (Intelligent-Tiering / InvalidObjectState) ────────
def object_tier(key: str, s3=None):
    """Return (storage_class, is_restoring, restored_ready) for a key.

    A HEAD always works even on archived objects. Use this before GET to avoid the
    'InvalidObjectState: operation not valid for the object's access tier' error.
    """
    s3 = s3 or _client()
    h = s3.head_object(Bucket=BUCKET, Key=key)
    sc = h.get("StorageClass", "STANDARD")
    restore = h.get("Restore", "")  # e.g. 'ongoing-request="true"' or '...="false", expiry-date=...'
    is_restoring = 'ongoing-request="true"' in restore
    restored_ready = 'ongoing-request="false"' in restore
    return sc, is_restoring, restored_ready


def ensure_retrievable(key: str, s3=None, restore_days=3, tier="Standard"):
    """If key is in an archive tier and not yet restored, issue a restore and report.

    Intelligent-Tiering archived objects need a restore request before GET. Returns
    'ready' | 'restoring' | 'requested'. Poll with object_tier() until 'ready'.
    """
    s3 = s3 or _client()
    sc, restoring, ready = object_tier(key, s3=s3)
    archived_classes = {"GLACIER", "DEEP_ARCHIVE", "INTELLIGENT_TIERING"}
    if sc not in archived_classes or ready:
        return "ready"
    if restoring:
        return "restoring"
    try:
        if sc == "INTELLIGENT_TIERING":
            # IT Archive Access / Deep Archive Access: the restore request takes
            # NO Days and NO retrieval Tier (passing either is rejected). The object
            # returns to the Frequent Access tier PERMANENTLY (not a temporary N-day
            # copy like Glacier), and Intelligent-Tiering has no retrieval fees.
            s3.restore_object(Bucket=BUCKET, Key=key, RestoreRequest={})
        else:
            # S3 Glacier Flexible Retrieval / Deep Archive: temporary restore,
            # requires Days and a retrieval Tier (Standard / Bulk / Expedited).
            s3.restore_object(
                Bucket=BUCKET, Key=key,
                RestoreRequest={"Days": restore_days,
                                "GlacierJobParameters": {"Tier": tier}})
        return "requested"
    except ClientError as e:
        code = e.response["Error"]["Code"]
        if code == "RestoreAlreadyInProgress":
            return "restoring"
        raise


def is_ready(key: str, s3=None) -> bool:
    """True if key can be GET now (non-archived, or archived-and-restored)."""
    s3 = s3 or _client()
    sc, _restoring, ready = object_tier(key, s3=s3)
    return sc not in {"GLACIER", "DEEP_ARCHIVE", "INTELLIGENT_TIERING"} or ready


def batch_restore(keys, s3=None, **kw) -> dict:
    """Issue a restore for every key up front (returns {key: state}).

    Restore is asynchronous, so kick them ALL off first, then poll with
    wait_until_restored -- do not restore-then-wait one file at a time.
    """
    s3 = s3 or _client()
    status = {}
    for k in keys:
        try:
            status[k] = ensure_retrievable(k, s3=s3, **kw)
        except ClientError as e:
            status[k] = "error:" + e.response["Error"]["Code"]
    return status


def wait_until_restored(keys, s3=None, poll=300, timeout=6 * 3600) -> dict:
    """Poll until every key is retrievable (or timeout). Returns {ready, pending}.

    IT Archive Access typically restores within 3-5 h (Deep Archive Access ~12 h),
    so poll on the order of minutes, not seconds.
    """
    s3 = s3 or _client()
    pending = set(keys)
    t0 = time.time()
    while pending and time.time() - t0 < timeout:
        for k in list(pending):
            if is_ready(k, s3=s3):
                pending.discard(k)
        if pending:
            time.sleep(poll)
    return {"ready": [k for k in keys if k not in pending],
            "pending": sorted(pending)}


# ── read one peak track ──────────────────────────────────────────────────────
def read_peaks(key: str, autosomes_only=True, s3=None,
               wait_restore=False, poll=30, timeout=1800) -> pd.DataFrame:
    """GET one peak file into a DataFrame (only PEAK_USECOLS, autosomes by default).

    Handles the archived tier: if wait_restore, blocks until the object is restored
    (polling object_tier), else raises a clear message telling you to restore first.
    """
    s3 = s3 or _client()
    state = ensure_retrievable(key, s3=s3)
    if state != "ready":
        if not wait_restore:
            raise RuntimeError(
                f"{key} is in archive tier (state={state}); call ensure_retrievable() "
                f"and wait, or pass wait_restore=True.")
        t0 = time.time()
        while time.time() - t0 < timeout:
            _, _, ready = object_tier(key, s3=s3)
            if ready:
                break
            time.sleep(poll)
        else:
            raise TimeoutError(f"{key} not restored within {timeout}s")

    raw = s3.get_object(Bucket=BUCKET, Key=key)["Body"].read()
    with gzip.open(io.BytesIO(raw), "rt") as fh:
        df = pd.read_csv(fh, usecols=lambda c: c in PEAK_USECOLS)
    df = df.rename(columns={"seqnames": "chrom", "Assay title": "assay",
                            "Target of assay": "target", "Biosample_name": "tissue"})
    if autosomes_only:
        auto = {f"chr{i}" for i in range(1, 23)}
        df = df[df["chrom"].isin(auto)]
    return df


if __name__ == "__main__":
    # smoke test on an S3-access instance: build the chip/atac/DNase catalog and
    # read the first non-archived track, without scanning every file.
    s3 = _client()
    cat = build_catalog(("chip", "atac", "DNase"))
    cat.to_parquet("atlas_catalog.parquet", index=False)
    print(cat.head())
    if "s3_key" in cat.columns and len(cat):
        k = cat["s3_key"].iloc[0]
        sc, restoring, ready = object_tier(k, s3=s3)
        print(f"first track {k}\n  tier={sc} restoring={restoring} ready={ready}")
