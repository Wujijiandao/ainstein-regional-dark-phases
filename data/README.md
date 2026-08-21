# SPARC RAR data provenance

The manuscript uses the public SPARC “Radial Acceleration Relation: All Data” table maintained by the SPARC collaboration. The SPARC webpage makes the table downloadable for scientific use and asks users to cite the SPARC master paper and relevant RAR papers.

This repository treats the table as **third-party data** rather than MIT-licensed code. For a public GitHub repository, the preferred packaging is to omit `RAR.mrt` and run:

```bash
python scripts/fetch_sparc.py
```

The manuscript snapshot has SHA-256:

`24aa7059dab7fa44787f7c11191052489899819370f6508621674769f3b72833`

If upstream content changes, inspect the change before updating the published regression values.
