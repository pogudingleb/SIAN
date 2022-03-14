array=(Goodwin HIV2 QWWC SEIAJRC) #NewBoneSimplified NewBrainReduced NewPCSK9 NewSuppl NewBoneFull)
for example in "${array[@]}"; do
    echo $example
    mkdir -p new_logs/$example
    /local/maple2021/bin/maple $example.mpl >new_logs/$example.out
done
# NewAtopicDematitis
