## Compare different-masked A.tha reference genome.

### Masked region size.

| Name                                                   | chrLength | masked size | coverage |
|:-------------------------------------------------------|:----------|:------------|:---------|
| Current version of soft_masked by ENSEMBL (version 57) | 119146348 | 38123558    | 0.3200   |
| Old version of soft_masked by ENSEMBL (version 52)     | 119146348 | 28109149    | 0.2359   |
| RepeatMasked with CONS-Dfam_withRBRM                   | 119146348 | 21267278    | 0.1785   |
| RepeatMasked with CONS-Dfam (version 3.7)              | 119146348 | 2393023     | 0.0201   |
| Not masked                                             | 119146348 | 186207      | 0.0016   |

RepeatMasked with CONS-Dfam_withRBRM:

```shell
#==================================================
#               number of      length   percentage
#               elements*    occupied  of sequence
#--------------------------------------------------
#Retroelements         8876      9977785 bp    8.34 %
#   SINEs:              583        87954 bp    0.07 %
#   Penelope:             0            0 bp    0.00 %
#   LINEs:             2346      1433309 bp    1.20 %
#    CRE/SLACS            0            0 bp    0.00 %
#     L2/CR1/Rex          0            0 bp    0.00 %
#     R1/LOA/Jockey       0            0 bp    0.00 %
#     R2/R4/NeSL          0            0 bp    0.00 %
#     RTE/Bov-B           0            0 bp    0.00 %
#     L1/CIN4          2276      1412082 bp    1.18 %
#   LTR elements:      5947      8456522 bp    7.07 %
#     BEL/Pao             0            0 bp    0.00 %
#     Ty1/Copia        1982      1902097 bp    1.59 %
#     Gypsy/DIRS1      3277      6474086 bp    5.41 %
#       Retroviral        0            0 bp    0.00 %
#
#DNA transposons      10551      5983962 bp    5.00 %
#   hobo-Activator     1540       567872 bp    0.47 %
#   Tc1-IS630-Pogo      915       237966 bp    0.20 %
#   En-Spm                0            0 bp    0.00 %
#   MULE-MuDR          4419      3313150 bp    2.77 %
#   PiggyBac              0            0 bp    0.00 %
#   Tourist/Harbinger  1002       336925 bp    0.28 %
#   Other (Mirage,        0            0 bp    0.00 %
#    P-element, Transib)
#
#Rolling-circles       3787      2282692 bp    1.91 %
#
#Unclassified:          269        77280 bp    0.06 %
#
#Total interspersed repeats:    16039027 bp   13.40 %
#
#
#Small RNA:            1167       143539 bp    0.12 %
#
#Satellites:           3260       978658 bp    0.82 %
#Simple repeats:      34152      1333320 bp    1.11 %
#Low complexity:       8606       414179 bp    0.35 %
#==================================================
```

RepeatMasked with CONS-Dfam (version 3.7):

```shell
#==================================================
#               number of      length   percentage
#               elements*    occupied  of sequence
#--------------------------------------------------
#Retroelements            0            0 bp    0.00 %
#   SINEs:                0            0 bp    0.00 %
#   Penelope:             0            0 bp    0.00 %
#   LINEs:                0            0 bp    0.00 %
#    CRE/SLACS            0            0 bp    0.00 %
#     L2/CR1/Rex          0            0 bp    0.00 %
#     R1/LOA/Jockey       0            0 bp    0.00 %
#     R2/R4/NeSL          0            0 bp    0.00 %
#     RTE/Bov-B           0            0 bp    0.00 %
#     L1/CIN4             0            0 bp    0.00 %
#   LTR elements:         0            0 bp    0.00 %
#     BEL/Pao             0            0 bp    0.00 %
#     Ty1/Copia           0            0 bp    0.00 %
#     Gypsy/DIRS1         0            0 bp    0.00 %
#       Retroviral        0            0 bp    0.00 %
#
#DNA transposons          0            0 bp    0.00 %
#   hobo-Activator        0            0 bp    0.00 %
#   Tc1-IS630-Pogo        0            0 bp    0.00 %
#   En-Spm                0            0 bp    0.00 %
#   MULE-MuDR             0            0 bp    0.00 %
#   PiggyBac              0            0 bp    0.00 %
#   Tourist/Harbinger     0            0 bp    0.00 %
#   Other (Mirage,        0            0 bp    0.00 %
#    P-element, Transib)
#
#Rolling-circles          0            0 bp    0.00 %
#
#Unclassified:            0            0 bp    0.00 %
#
#Total interspersed repeats:           0 bp    0.00 %
#
#
#Small RNA:               0            0 bp    0.00 %
#
#Satellites:              0            0 bp    0.00 %
#Simple repeats:      40167      1627440 bp    1.36 %
#Low complexity:      10698       583930 bp    0.49 %
#==================================================
```

### Segmental duplications.

| Name                                                   | chrLength | defined size | coverage | copy2 count | BISER  | BISER merged | ASGART |
|:-------------------------------------------------------|:----------|:-------------|:---------|:------------|:-------|:-------------|:-------|
| Current version of soft_masked by ENSEMBL (version 57) | 119146348 | 9776236      | 0.0818   | 2625        | 2561   | 1017         | 238    |
| Old version of soft_masked by ENSEMBL (version 52)     | 119146348 | 11434422     | 0.0957   | 3155        | 4094   | 1580         | 418    |
| RepeatMasked with CONS-Dfam_withRBRM                   | 119146348 | 13098118     | 0.1096   | 3480        | 6339   | 2036         | 670    |
| RepeatMasked with CONS-Dfam (version 3.7)              | 119146348 | 19236669     | 0.1610   | 5019        | 95566  | 7247         | 9341   |
| Not masked                                             | 119146348 | 19188917     | 0.1606   | 4987        | 100826 | 7417         | 11775  |

