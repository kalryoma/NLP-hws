for i in {1..10}
do
    ld=$(python -c "print ${i}*1.0/10")
    score=$(python answer/mybaseline.py -l ${ld}|python score-segments.py)
    echo $ld
    echo $score
done
