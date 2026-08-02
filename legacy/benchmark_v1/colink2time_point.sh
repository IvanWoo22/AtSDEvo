METHOD=$1
REF=$2

cut -f 4,5,6,7,8 "${METHOD}"_"${REF}"/time-point.tsv | sort | uniq -c \
	| awk -va=1 '{
				count0 = 0;
				count1 = 0;
				count2 = 0;
				count3 = 0;
				for (METHOD = (a+1); METHOD <= NF; METHOD++) {
					if ($METHOD == 0) {
							count0++;
					} else if ($METHOD == 1) {
							count1++;
					} else if ($METHOD == 2) {
							count2++;
					} else if ($METHOD == 3) {
							count3++;
					}
				}
				if (count0 >= 0 && count1 == 0 && count2 == 0 && count3 == 0) {
					print;
				}
			}' | awk '{print $2 $3 $4 $5 $6 "\t" 0 "\t" $1}' >time1-0.tmp
cut -f 4,5,6,7,8 "${METHOD}"_"${REF}"/time-point.tsv | sort | uniq -c \
	| awk -va=1 '{
				count0 = 0;
				count1 = 0;
				count2 = 0;
				count3 = 0;
				for (METHOD = (a+1); METHOD <= NF; METHOD++) {
					if ($METHOD == 0) {
							count0++;
					} else if ($METHOD == 1) {
							count1++;
					} else if ($METHOD == 2) {
							count2++;
					} else if ($METHOD == 3) {
							count3++;
					}
				}
				if (count0 >= 0 && count1 > 0 && count2 == 0 && count3 == 0) {
					print;
				}
			}' | awk '{print $2 $3 $4 $5 $6 "\t" 1 "\t" $1}' >time1-1.tmp
cut -f 4,5,6,7,8 "${METHOD}"_"${REF}"/time-point.tsv | sort | uniq -c \
	| awk -va=1 '{
				count0 = 0;
				count1 = 0;
				count2 = 0;
				count3 = 0;
				for (METHOD = (a+1); METHOD <= NF; METHOD++) {
					if ($METHOD == 0) {
							count0++;
					} else if ($METHOD == 1) {
							count1++;
					} else if ($METHOD == 2) {
							count2++;
					} else if ($METHOD == 3) {
							count3++;
					}
				}
				if (count0 >= 0 && count1 == 0 && count2 > 0 && count3 == 0) {
					print;
				}
			}' | awk '{print $2 $3 $4 $5 $6 "\t" 2 "\t" $1}' >time1-2.tmp
cat time1-{0,1,2}.tmp >"${METHOD}"_"${REF}"/time1.txt
for TIME in {2..5}; do
	cut -f 4,5,6,7,8 "${METHOD}"_"${REF}"/time-point.tsv | sort | uniq -c \
		| awk -va="$TIME" '{
					count0 = 0;
					count1 = 0;
					count2 = 0;
					count3 = 0;
					for (TIME = (a+1); TIME <= NF; TIME++) {
						if ($TIME == 0) {
								count0++;
						} else if ($TIME == 1) {
								count1++;
						} else if ($TIME == 2) {
								count2++;
						} else if ($TIME == 3) {
								count3++;
						}
					}
					if (count0 >= 0 && count1 == 0 && count2 == 0 && count3 == 0) {
						print;
					}
				}' | awk -va="$TIME" '$a==3{print $2 $3 $4 $5 $6 "\t" 0 "\t" $1}' >time"${TIME}"-0.tmp
	cut -f 4,5,6,7,8 "${METHOD}"_"${REF}"/time-point.tsv | sort | uniq -c \
		| awk -va="$TIME" '{
					count0 = 0;
					count1 = 0;
					count2 = 0;
					count3 = 0;
					for (TIME = (a+1); TIME <= NF; TIME++) {
						if ($TIME == 0) {
								count0++;
						} else if ($TIME == 1) {
								count1++;
						} else if ($TIME == 2) {
								count2++;
						} else if ($TIME == 3) {
								count3++;
						}
					}
					if (count0 >= 0 && count1 > 0 && count2 == 0 && count3 == 0) {
						print;
					}
				}' | awk -va="$TIME" '$a==3{print $2 $3 $4 $5 $6 "\t" 1 "\t" $1}' >time"${TIME}"-1.tmp
	cut -f 4,5,6,7,8 "${METHOD}"_"${REF}"/time-point.tsv | sort | uniq -c \
		| awk -va="$TIME" '{
					count0 = 0;
					count1 = 0;
					count2 = 0;
					count3 = 0;
					for (TIME = (a+1); TIME <= NF; TIME++) {
						if ($TIME == 0) {
								count0++;
						} else if ($TIME == 1) {
								count1++;
						} else if ($TIME == 2) {
								count2++;
						} else if ($TIME == 3) {
								count3++;
						}
					}
					if (count0 >= 0 && count1 == 0 && count2 > 0 && count3 == 0) {
						print;
					}
				}' | awk -va="$TIME" '$a==3{print $2 $3 $4 $5 $6 "\t" 2 "\t" $1}' >time"${TIME}"-2.tmp
	cat time"${TIME}"-{0,1,2}.tmp >"${METHOD}"_"${REF}"/time"${TIME}".txt
done
cut -f 4,5,6,7,8 "${METHOD}"_"${REF}"/time-point.tsv | sort | uniq -c \
	| awk '$NF==3{print $2 $3 $4 $5 $6 $7 "\t" 0 "\t" $1}' >"${METHOD}"_"${REF}"/time6.txt
rm time*.tmp
awk -va="${METHOD}" -vb="${REF}" '{print b "\t" a "\t" ($3-$2) "\n" b "\t" a "\t" ($6-$5)}' \
	"${METHOD}"_"${REF}"/links.2copy.tsv >>links2_length_distribution.tsv
for TIME in {1..6}; do
	awk -va="${METHOD}" -vb="${REF}" -vc="${TIME}" '{sum+=$3};END{print b "\t" a "\t" c "\t" sum}' \
		"${METHOD}"_"${REF}"/time"${TIME}".txt >>links2_timedistribution.tsv
done
awk -F"\t" -va="${METHOD}" -vb="${REF}" '$1==b&&$2==a{sum+=$4};END{print b "\t" a "\t" sum}' \
	links2_timedistribution.tsv >>links2_timedistributionallcount.tsv
wc -l <"${METHOD}"_"${REF}"/links.2copy.tsv >>links2_allcount.tsv
