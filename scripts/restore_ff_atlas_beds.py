#!/usr/bin/env python
"""DEPRECATED -- do not use. Superseded by scripts/fetch_ff_atlas_beds.py.

This script was written on a false premise: that the 73-tissue FF openness atlas
was built from ENCODE objects whose accessions survived only in an archived-S3
metadata index, and therefore that regenerating the source BEDs required AWS
credentials and a Glacier-tier restore.

That premise is WRONG. The 73-tissue atlas was built ENTIRELY from the PUBLIC
ENCODE portal (an hg19 DNase-seq narrowPeak query), curated to classification==
"tissue" plus a fixed ADDITIONS list of immune / myeloid-progenitor / cfDNA-
relevant biosamples. Every source BED is re-fetchable from
https://www.encodeproject.org with no S3 bucket, no IAM role, and no restore --
the accessions are recovered by re-running the same public query, not read from
any private index. (The 73rd track, neutrophil, is the one exception: ENCODE has
no hg19 neutrophil DNase-seq, so it comes from local BEDs -- see the fetch
script's note.)

USE INSTEAD:

    python scripts/fetch_ff_atlas_beds.py --dry-run   # preview the 72-biosample panel
    python scripts/fetch_ff_atlas_beds.py             # download BEDs + manifest -> data/access_ff/

This file is kept only as a tombstone so the corrected provenance is on record.
"""

import sys

if __name__ == "__main__":
    sys.exit(
        "restore_ff_atlas_beds.py is DEPRECATED (its S3-restore premise was "
        "incorrect).\nUse scripts/fetch_ff_atlas_beds.py -- the atlas BEDs are "
        "public ENCODE, no S3/AWS required."
    )
