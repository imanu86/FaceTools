# Status and roadmap

Prototype. Native tests cover the geometry core; GUI acceptance is partial (see TEST_CHECKLIST.md). Known open: "Cut and extend" needs verification on U/Y/T shapes, ambiguous guides, invalidation and Undo/Redo; strict-export modal; project save/reopen. Out of scope: shortening, general CAD Move Face, curved surfaces/BRep, free points across arbitrary instances, notarised distribution, CI. Undo/BBS does not serialise `exterior_facets` (conservative global block kept).
