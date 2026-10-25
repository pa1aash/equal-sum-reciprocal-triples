# Data

Every data file of this repository: 42 files, 11.42 MB. For each, its size, its SHA-256 checksum, the script that produces it and the element of the paper it supports. `make data-check` compares every file with the checksums in `SHA256SUMS`. Every step of `make` also fails if it changes one of these files: regenerating committed data is a reproduction test, and the result must be byte-identical.

Licence: CC-BY-4.0 (`LICENSE-DATA`).

| file | size | SHA-256 | produced by | paper element |
|---|---|---|---|---|
| `checks/classes4800/check_descent.txt` | 2.8 KB | `d2a329c963fcbb5bdccd9ae73fbbd57ee221bded1d31c564bac51484bbe77ea5` | `checks/classes4800/check_descent.py` | the 2-isogeny descents of Section 5 |
| `checks/classes4800/check_enum.txt` | 2.2 KB | `c2006ef3d8fef4eaaafd212eb9eda9997b0d754e009d6754b7dc24e61f52ee58` | `checks/classes4800/check_enum.py` | Section 3 (orders of $P'\pm P$), Section 6.1 (histogram, admissible triples) |
| `checks/classes4800/check_groups.txt` | 5.8 KB | `25aef4f6a8ba53d9bf779b718ddbacfee96aa05cbc057ebafed078ffa1959a96` | `checks/classes4800/check_groups.py` | Theorem 4 ($nP\ne O$), Lemma 2, Theorem 1, Proposition 2 |
| `checks/classes4800/check_isosceles.txt` | 1012 B | `54a6b394536875e6a321aa176b6206ae9631bd98efc074ace3684ed29f4d5078` | `checks/classes4800/check_isosceles.py` | Theorem 2 (the copies with $k\ge4$, a subset of those counted in Section 4) |
| `checks/classes4800/check_theorem5_family.txt` | 1.7 KB | `ba3fa3160b5649b07202fe7c8d9fb74b3d027a553008e45baf930237fab348b1` | `checks/classes4800/check_theorem5_family.py` | Lemma 4 and Theorem 3 |
| `checks/classes4800/enum_4800.txt` | 4.58 MB | `ced0f705471a2f6efd49e0260e9bd0bbccd95a9ab45783be24623d0d6fd0a52c` | `checks/classes4800/enum_classes.c` (`enum_classes 10 4800 4800`) | the second enumeration of Section 6.1 |
| `checks/enumeration6000/check_enum.txt` | 5.3 KB | `45bd0c9ccdda6ecf80db001741f8d40ca8cb2eae0298034b4c2e38c969def814` | `checks/enumeration6000/check_enum.py` (standard output) | Table 2 (counts and local exponents to $S=6000$), Section 6.3 |
| `checks/enumeration6000/data/fibres_10_2000.txt` | 958.5 KB | `e3a54cece65650fccb6299c364dcf36509a44bac39d4522e9b09b5c554b5e051` | `checks/enumeration6000/enum.c` (via `check_enum.py`) | class listing to $S=2000$, compared with the pure-Python enumeration |
| `checks/enumeration6000/data/summary_10_6000.txt` | 172.7 KB | `b9d7c5f70b31a8ff0d473a18211444e253ab202637c63743f813c0909e3b2cf5` | `checks/enumeration6000/enum.c` (via `check_enum.py`) | Table 2, Section 6.3 |
| `checks/enumeration6000/tab_enum_reference.txt` | 1.4 KB | `603dcc8bcd3981d2e3ba5aa14e7387e84443e90edb3bb5a48795b15fcc5a2ac2` | `checks/enumeration6000/check_enum.py` | listing of the classes of small sum (by-product) |
| `checks/note/check_note.txt` | 6.0 KB | `bc37f163b6851885b6d6e76f21b7c090ed47c28b937699c5716cde648bd81505` | `checks/note/check_note.py --enum` | Lemmas 1 and 3, Proposition 2, Theorems 1 and 2, Proposition 4, Table 1, Table 2 (row 6000) |
| `checks/point3P/check_point3P.txt` | 517 B | `5af28a657293d4b4337d180c0d947b769633bd6a08a7b4e36a272c4a6568d419` | `checks/point3P/check_point3P.py` | the point $3P$ of Section 5 |
| `data/degeneracies.csv` | 16.7 KB | `481f28ea7949d3a5090e907cad2568cd1b5606a00ed94daeb38460a39a95eb55` | `code/enumerate_degeneracies.py` | reference counts for $S\le600$, read by `check_committed.py` |
| `data/degeneracy-groups.csv` | 173.7 KB | `ee97813e67c7f74997066125f9ab839e1e4ba7d7d702402ba184f47d186c48a8` | `code/enumerate_degeneracies.py` | reference classes for $S\le600$, checked by `enumerate_fast.py` |
| `figures/out/F9.pdf` | 59.4 KB | `dff8c0428cd73f16c07a6ebdc34255ea8ea6d4077a2b9457b7da7d116b980e9d` | `figures/src/F9.py` | Figure 1 |
| `figures/out/F9.png` | 36.1 KB | `0d8015b244f2a6068f0773b53f5da6c305ee98b37f0ec99d0ba82995851fa69c` | `figures/src/F9.py` | Figure 1 (preview) |
| `figures/style/palette.json` | 2.1 KB | `5be6369e62e68269866b7cf419671b39e54b9da21376bff53f369abb529225b6` | fixed input (colour palette) | Figure 1 (style) |
| `figures/style/third_party/SOURCES.json` | 749 B | `1f481b8e1aea83df02a8779d2ffe5d30aeb7b4218799a3c32087a54f5f9f4c5e` | fetched colour tables (`SOURCES.json` gives their sources and licences) | Figure 1 (colour maps) |
| `figures/style/third_party/crameri_acton.txt` | 6.8 KB | `068da7de02af4d0607faa08b233c26065a9c582911e437bc686962ebe021ceda` | fetched colour tables (`SOURCES.json` gives their sources and licences) | Figure 1 (colour maps) |
| `figures/style/third_party/crameri_batlow.txt` | 6.8 KB | `918dd393177e126685437f234f4c1f7c8a446437c1c68820f3bc1a8423e91bd1` | fetched colour tables (`SOURCES.json` gives their sources and licences) | Figure 1 (colour maps) |
| `figures/style/third_party/crameri_bilbao.txt` | 6.8 KB | `b641316b9bc2a952abecdf47494ae857d85017fb3715dfba68e6836fc8c6946f` | fetched colour tables (`SOURCES.json` gives their sources and licences) | Figure 1 (colour maps) |
| `figures/style/third_party/crameri_davos.txt` | 6.8 KB | `0ca57a52e94c12c05aef0abf22f363f9fed0cb831cdde2aed76ee8bb4d950bd9` | fetched colour tables (`SOURCES.json` gives their sources and licences) | Figure 1 (colour maps) |
| `figures/style/third_party/crameri_devon.txt` | 6.8 KB | `b9573d181e49ead6e6daef40641b82b04116e2f9dacb94b4a9fdd3afc2f55dd0` | fetched colour tables (`SOURCES.json` gives their sources and licences) | Figure 1 (colour maps) |
| `figures/style/third_party/crameri_glasgow.txt` | 6.8 KB | `1fd603c787431927b006a99c8590cdaa77a285da004f017dc766b943332246ce` | fetched colour tables (`SOURCES.json` gives their sources and licences) | Figure 1 (colour maps) |
| `figures/style/third_party/crameri_imola.txt` | 6.8 KB | `aa0689d10ae651cf0d04c23fddf1788e7a3922a76c415414e92d3f9b160926d5` | fetched colour tables (`SOURCES.json` gives their sources and licences) | Figure 1 (colour maps) |
| `figures/style/third_party/crameri_lajolla.txt` | 6.8 KB | `7541ff8f5267045b8fc716c53d2fa9153882fe23471c74f8c3692ff57aef1971` | fetched colour tables (`SOURCES.json` gives their sources and licences) | Figure 1 (colour maps) |
| `figures/style/third_party/crameri_lipari.txt` | 6.8 KB | `59a26ca79c66da1674f21e4ef83dd6a428c7f1e768327fc4b09e3242b7b81829` | fetched colour tables (`SOURCES.json` gives their sources and licences) | Figure 1 (colour maps) |
| `figures/style/third_party/crameri_oslo.txt` | 6.8 KB | `4e29e907d1a07e05481ee6856efe4b5266d85f89eea8c82affde764aacaf1f89` | fetched colour tables (`SOURCES.json` gives their sources and licences) | Figure 1 (colour maps) |
| `figures/style/third_party/machado2009_cvd.json` | 658 B | `f8307014f929ebd7c1f61a2340940dcf8fe7a8e9a2db66916cd5173fba36e554` | fetched colour tables (`SOURCES.json` gives their sources and licences) | Figure 1 (colour maps) |
| `figures/style/third_party/viridis.txt` | 8.2 KB | `01aa44c85dc376fadfc8595072aed0edffb155b4edb1090381064fc63d0b73bc` | fetched colour tables (`SOURCES.json` gives their sources and licences) | Figure 1 (colour maps) |
| `tables/density.tex` | 686 B | `f92ae7889d80156727afb0365dba3fefadde07dc29046c8ec6ea6bad0d4e07a6` | `tools/make_tables.py` | Tables 1 and 2 |
| `tables/fibres.tex` | 689 B | `f17a7999d7a41377632eab970fa51c4955858f36700fad09009522f0c42ab846` | `tools/make_tables.py` | Tables 1 and 2 |
| `theory/diophantine/data/birthday.csv` | 1.3 KB | `8d479ba1b5c374f533d2c067437825ed5d710c1167cd3dd0549ec5fda1f385a1` | `theory/diophantine/growth_fits.py` | input of the model comparison of Section 6.3 |
| `theory/diophantine/data/exponent_fits.txt` | 4.3 KB | `a32b9080400a0a0874660a895b00d278c7f7fb33f8a24fae54eeb330a3d8eed6` | `theory/diophantine/exponent_fits.py` | the fits of Section 6.3, Figure 1(b) |
| `theory/diophantine/data/families.txt` | 4.2 KB | `e920e247a6b41a2218e0a48e04b5c8e677936274c8eb8daa03d1bfaeec0fda24` | `theory/diophantine/families.py` | dual pairs (Section 3, 6.1), isosceles count (Theorem 2), lines (Section 6.4), Figure 1 |
| `theory/diophantine/data/groups.csv` | 5.02 MB | `d56e8875af1bbb299339a8728435a3ee33bd0ce23440f220f475011d88eaebb1` | `theory/diophantine/enumerate_fast.py` | every class with $S\le4800$: Table 1, Theorem 4, Section 6, Figure 1 |
| `theory/diophantine/data/growth_fits.txt` | 4.6 KB | `0965cbe9933a732de16d33afd376e9aa3d5a3694cba162e661087dc16ce40325` | `theory/diophantine/growth_fits.py` | the model comparison of Section 6.3 |
| `theory/diophantine/data/line_families.csv` | 1.1 KB | `a610e5aae4ba40007f7cc8931c581bc919634059cf3d4655f3b0ad354b68bada` | `theory/diophantine/families.py` | the 189 lines of Section 6.4 |
| `theory/diophantine/data/per_S.csv` | 319.3 KB | `b9da0d34b414186ebaf84bc4923e0cebde357c75fbe9b8dd43dd9a33e04dea9f` | `theory/diophantine/enumerate_fast.py` | counts per sum: Table 2, Sections 3 and 6, Figure 1 |
| `theory/diophantine/data/ranks.txt` | 2.6 KB | `b10991112c38cea64ea437ad62477502c3336a50714e73de12cd82d74b235fe1` | `theory/diophantine/ranks.py` (PARI) | ranks and torsion of Table 1 and Section 5 |
| `theory/diophantine/data/torsion_pairs.txt` | 548 B | `99948aaef877239ff43158590948f4e9fdb3d2fee27afb19709c44b326ec5a2f` | `theory/diophantine/torsion_pairs.py` | Section 3 (pairs that differ by a torsion point) |
| `theory/diophantine/data/variety_checks.txt` | 1.7 KB | `d6ae677a8ba21f4edde7bfce2dd9121c945de9c1709ecd64f3e03efabd6e75fd` | `theory/diophantine/variety_checks.py` | Propositions 1-3 (pencil, maps, base points, reciprocation) |
