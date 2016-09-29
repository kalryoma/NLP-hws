for i in {10..20}
do
    ld=$(python -c "print ${i}*1.0/100")
    score=$(python answer/mybaseline.py -l ${ld}|python score-segments.py)
    echo $ld
    echo $score
done
