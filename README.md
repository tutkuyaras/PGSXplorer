# PGSExplorer
PGSExplorer is a bioinformatics workflow designed to calculate polygenic scores by processing genomic data through quality control steps. Optionally, it can utilize tools such as [PLINK](https://www.cog-genomics.org/plink/), [PRSice-2](https://choishingwan.github.io/PRSice/), [LD-Pred2 (grid)](https://privefl.github.io/bigsnpr/articles/LDpred2.html), [LD-Pred (auto)](https://privefl.github.io/bigsnpr/articles/LDpred2.html), [PRS-CSx](https://github.com/getian107/PRScsx), and [MUSSEL](https://github.com/Jin93/MUSSEL). The workflow requires genomic files in PLINK format (.bed, .bim, .fam) and GWAS summary statistics for two different populations as input to complete the analysis.

![](https://viewer.diagrams.net/tags=%7B%7D&lightbox=1&target=blank&highlight=0000ff&edit=_blank&layers=1&nav=1&title=Ba%C5%9Fl%C4%B1ks%C4%B1z%20Diyagram.drawio#R%3Cmxfile%3E%3Cdiagram%20name%3D%22Sayfa%201%22%20id%3D%22A8VOpgxVlnepOkoelsrt%22%3E7V1tk5s4Ev41U7f7wRRCvH60xzPJ1WW35jK7l7v7ksIg22wweAHPy%2F76lUACAcLGNtg4YVIVGyEL6H70dKvVEnfwfvP2IbK3619CF%2Fl3iuy%2B3cH5naIAQ9HxByl5z0p0y8wKVpHn0kpFwbP3F6KFMi3deS6KSxWTMPQTb1sudMIgQE5SKrOjKHwtV1uGfvmqW3uFagXPju3XS794brLOSk3FKMo%2FIm%2B1ZlcGupWdWdjOt1UU7gJ6vSAMUHZmY7Nm6DPGa9sNX7ki%2BHAH76MwTLJvm7d75BOxMollv3tsOJvfcoSCpM0P%2FqNMF2v45%2FMf%2F58Ev1uTrfdbsp4AFWbtvNj%2BjgqD3m7yzqRDHnBbvxy9gxcUJehNpCZ7wVoonhdDCIUblETvuB79lS5TEVH0qLKWHb8WutBplTWnBp0hgGp%2FlbdciAF%2FoZI4SirqYamQ5%2FYwhD7ZC%2BQ%2FhbGXeGGATy3CJAk3d3C2Tjb4qnOAv7K6U99bkTpJuMWlNj1ysEhRhAviJAq%2FofvQD6P0GhDq0IL4TmdLz%2Fe5cuNhqj%2FI5Bdre0vuZ%2FO2Ij1Tsv%2FaRUhy%2FHDnfo1R9OI5KP6K1bH0VrvIJnf4FTeFfzHbhh656sMLvnhMb7Ou4wOgqaq%2BUcVAK2l4AoyahqFAw6ysew1r5vVxj%2B%2BiDHwA1ZpYVEY1JblApS%2FBAOtkwRyQ9NnyAhVpyUAgLQGKFNgbUShax9LKxX%2B6tHhxKEAy4CUlYskdCySXcCcCgXqDQOS%2BmEbXuxVILuDT5cHMrCVZ%2FJ%2BAkhukhXlLkvk%2F2JPwcq7lhaf7CbGIIRYTL0X9z13ITkzi1MWc4gpA2b6lpo2dx99W5PPR87H984IVcda8OM6%2BPf%2F6FLNL4FvOrpL9oKY2LP%2BkZrVTD1Bg23kLTuvscRA2nuumZvp17SXoeWs75Jqv2MbjstTfRERuMn3eR3vj%2BUStH%2BzI3oSBS8t5VyL9a2%2Fkc9y2tvFqmZ4NKKkCjHDw0gRcDSypL9Om1C2bGANe4Hovnruz%2FRwKuKCAQVa0qBVE1RJ8m7VauKzU2IipfZgyYcXkVzAl8pcuCir99PHTcSI6mt4VU9J5WSn1DtjI79htUHsSmNDxzrpGvLUDIaU7GeQInUerhf0TfnB8cbny8XMKRDnl%2FyUFL%2FkJxm4YZ2jPzxf2QabmQSY9b0K72DS9DTtKhJbjTtGqsRANS4CUpsP5%2FIhJREtlgkvm5DvRm0aeWsOCPFQX5HUZdk5qRimayQSfnylOZHLPT4jHt%2BQ0HeHirxWOISfJ0JacSYmIFID0MBuokmM2VF3YsedIWzvBXBZgIvY%2Fk8gOqUto7Tm%2FGdezV5%2FRS9ZKgtLW5xpXkTBlFrUhZ2TS1UvnMsJizwXd9C%2BvUzkLDKADPT3LiSzrjG2Ukp%2FK9VH0TI3QWl5TK%2FQnv3OaUixJ1spuGav2mj9ldi2uBdp9i2a4c7gP5%2BUr7m6q2EkPcwDxhWVY03o1%2FOcGLuvIDRYuj4MRQ%2BPa8Tq3OpztqkdCgAo1tR4J0VVDNw1c7lfCMLmRPBinaRmcCXeJ7wX4nlgIUhR%2FeUHI3kiundgszlIxmcv07wiT2T7UciCWosqS2Rentxj0nagGHhMVWYsJhKePgjwodRDiqNNGRhp1ykgVRumiBjxKFH0o8q2sRGrPLUnlacGsWXMgCn30NqwVjczO8XuaxHN0lEglPg3v%2BKhlFxJK1gEXEirigQlQOPn35R%2BpLeIF5U7RwKY1AmscD9SCyOmMipJ9kmbfWftExIB8C8LEWbPBQoWRXdPQVSGHTzUVgj2cuAlfqJZJrQhhL40%2FDhM74Y6x4UH8MXI9%2FtAPnW%2B5MMoYOiLWrTdFe%2FfzrwhBfQVJVOPGEeNqyHSFiDGVBdT1m0cGpZ76rNZlcdJiymPQOFkuHcey9s6T9aY7ZlyvprwW0zJxyacB%2B72npR%2B%2BOms8xJXwuDeyneRrGH3dIDveRWiTSq8hTpTY1HOz5F4EzibHDLEZ5xSgKOLoRV9zG1pzNHqtCgMXZMAwYfc2zWIVWkNMoSGqSOPcy0wLxQ8%2FIv8FkV5Xbws0BrD3hiXJMzSM2bBwvG3cFDvkUGbH22xwtPTeCF%2FMyphUqrHF3%2Bx1uLGz2I%2FjBasZHQ3M1aLsM1WtIur56pz8E5FSb8gEannexNRLDuYEShqPW60%2BryIijt48dU25DdQqI2r7RK0GJQBvCLWiuPoAUQtH1PaJWquCWmvYoBWljQ0QtOoI2l4dBIhRqw4bqaK8pQEiVRuR2itSjeEjVRQPHSBS9RGpfSKVxO2HjlRRHHaASDVGpPaKVG3wSLVEfuo4X3um4qFcG52YBqd1naHiWnO2bNGAgKC4lEQ%2BJ3JBHu%2BX6eOYvHjKBD5gM9gUENCUTJ4G6ssXLpwRyzIMOEB8Rr6dIDdAcYxP3K%2BR821U8wHCryxfg4pUClnUMzVEau5tRY9a7%2FWfEZlPJb8Ml%2Fi%2F%2BW7rYyUk6I7lv48K36two2zecb82eIVr9Yle45IaN1ukkDDj7G3SpbGHFyse9AfSBLxZvh5WlN%2BRXmzKfEBZ5BDS%2B5mvk4Ss%2FZ0SSSiPjhuokpeuY8SAiCSHpJ8%2BkjQ7%2FEHKY%2FIZBuyrBhR6YpJONU%2Fc8DWYLCZAMaUttmz0OjN%2B9S7Faki8yoRo1pBrDuox6QEnmAtd0q3GYKcGyrnvsl5DGQaiIPu9KO0caUaL7PfvEGnYOGLle7Y%2FsbfbCZEhSRl%2BjMnUeIzvfbKwo0k6uT4hZRzsmriqDkcXLe1d6nRl58LIRbXVv7XRUwW9fWK1BE4zXzDAL8UQZqH2RXotsqC%2BQyjeOukpsiYZw%2Be5H9OiXpnnHoF1z%2BV7XYPnFLm8ccH1eQ7WTe7HLw%2B4IBu3j777AcOly1J5TeEE6BIQuF3XG5ZbdRVLCyzKdJmbtPA29NvS3ozqPpBFCExJKdsXwdYjok0jTKO3JZCwnoj7EWGNhH%2B9r4gZeP9BIy9Hrw4maz1kWVE1SwYkA5ilurEVsPWOjDstgACohqbKlgmNSzqobTJ4vz8f4kIOajUv2fUi%2FAhZQRzuSL3%2BcAj0qhOLeef6XqvZIuF%2FRNxNIg67JBIEw4McEK4LHOAkuDlOgveJT6hUwakc3I3iolPgRosp8O%2BQG688op%2Fqpoy7Xncj%2BsbdsPYx5wEn8aLDeaPFBnIjDrvGoW5MZ9bjVXFIfEagH4AitNggtGS%2Frd6GplZzXvC4O8%2B4O8%2B4O8%2B4O4%2FIw%2FxRd%2Bdp3p202TPWy8kjhiUJJrguu38PEEwq1Mc3gTvNlu3PHd%2BOMT2UQVGOER4rQOSWdmjfu7uGKBzPyiKSv%2Ba98G2JBUav8ES2wea0Y2iSVtNPaeBSbhKP9CMH0VYKNQgaLtt6qJXbweZshZJaO6k6cymco%2BHmpWlNo%2Bc8slszsfe27%2Bz8LGHs6X7aOJSuIgjrOO3nlfFw6%2B0k%2BH0%2FCMjSPcxTAWqzO21%2BBOxg00zC8Zt4y5Ksm0JH7kwsVn4QLpcx6gccghTBmvIuvWlzTpfFFpX5Dk71vThEMziG1izxM18P0MV%2BZz2%2BHqAY3lzn9QBKg%2FabB0NlbmxpEXt8O0CLUfklFCwjUzbNuoIXMoJIvyEFy5UtM1ismuX8CLZGuqjCYYsEn2v26GLro6so%2FBQ3t%2BDtq2m1zV56%2Fbi17DUOQ3FrrXLQ09I6cWpVs%2FpGjooq%2B%2FZqhUmfR74fgIWYqk7uv%2B%2F%2FQXJudHtDum6wiLflqSDRrFG6OOqMa3Z281yWUJo5dGiuS7C%2F%2FXCeJc1zOuoJhjvwUJum0k4ZeADIxjGn8sMFBhfCZdFnkTB685L%2F0jPk%2B%2F9IuaTRo%2FkbV23%2Bzh08ocjboFTVDWROBZfxXxsHcSD8Xn33BmDbJx0dpqg0BM1LM%2FrBOMXpRPK0RkGYvJNI980FLPJe1AlvaGwrwwHzhtVih2lOKS2mNIUsgwvJ9Aajk4C9nZKMrv5ASfJOj%2B1dEpI4c5Ssw1UY4GFASHx6Lh7OJwcRWNEfArWRwOQzCQxjJnrPmZAccC2Rw6Kp9Ii1dQbxsTcjZjywp6bCnMpjKFL4lsduKFKr%2BKqqVWG2thSpVRIIYLWhvilSuHtDRxT54cv0%2BQbZMeeKs9kRk6MKy8GJG6DKFtM3I1UOmCrVQVGlrpanwKBZiaS3pkqlwrnVNKreqbI5neWHpcrOZr5IBnLl5RpKJ%2FirbJ5UCU71x6OK3CJe2PdLgWHFT8HikOrRcSB8z%2B0%2BUZ8XSbWad%2BFr3ZFg07Ds0z9%2F%2FReu8NPTnTL77ee2fYpFwqNsF6nCbJCedSAFu5TlIlxNVYm%2FT6179WFGumzkUGMF2VIo3njRdOxPaEl0orJcbra9leDlJ4%2FpH0cGfvrbtv25wGz7yTZLUkoLrMvTq7pkybrgXRWgSEApgc6Q8to94K6DAK%2FZgLvmbfIbkFZk%2BG8574KWcg4HjwpFBIE98OOdsMaXvpd8KLWOqbyFmoE50C96AF15NYpwZX8pywlIKuxv5t5qfvvk%2BTz2%2Bdlz0EQZCey6BKYBybB00xgKhYnWfXZEYUdjbaSwE2CXj14kU2war09s%2BR50fbwp%2FdN88hQhdyS2KxObpUuDoTVFFoXnO6K15pdqjLTWOa2ZlmSYYKCspvTqrk3u47dBkpo1vTdmjz8EqQEMIU0MwOvwWvOLh87mteb3rgyT1xqWSd4GrwGoSroBRAZzCMzWwSReE7P98vvz88OnkdiuTGymLqkD8tc6mApp4rXmt%2FSMvNY9rxmWpA10FAqO2HHI2UX%2B%2ByyynW9kcvSQPMu0VF1OSvt6U4Z%2BPo3GtBhlEmmeb6ur7%2FjXQVQ3mmCTnwe2PAHVfJIO8%2BG%2Fx1RMtnRjoKmYClsKcHYqJrh0KmaPk%2BcSFjwKXpBzgzPoapepmNhr0M%2FDXb7DUvkXfWZ0t1j%2BeHM0wnrbUGhEr6TOVHv%2FyTSiX5pGOkgd%2BM5ycPIO1EEOjpxnZ10%2FXfGO7kXCVS92IYEPfwM%3D%3C%2Fdiagram%3E%3C%2Fmxfile%3E)
