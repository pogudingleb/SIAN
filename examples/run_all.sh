array=("Tumor.mpl" "QWWC.mpl" "Fujita.mpl" "Goodwin.mpl" "HIV2.mpl" "New_AKT.mpl" "OralGlucose.mpl" "SEIQRDC.mpl" "SEIRP.mpl" "SIRAQJ.mpl" "Treatment.mpl" "new_SEIR.mpl") #(HIV2 Goodwin) #SSAAIR  SIAR Treatment Tumor) #Goodwin HIV2 QWWC SEIAJRC) #NewBoneSimplified NewBrainReduced NewPCSK9 NewSuppl NewBoneFull)

for example in "${array[@]}"; do
    echo ${example::-4}
    rm new_logs/${example::-4}.out
    mkdir -p new_logs/${example::-4}
    /local/maple2021/bin/maple $example >>new_logs/${example::-4}.out
done
# NewAtopicDematitis
