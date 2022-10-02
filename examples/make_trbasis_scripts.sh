array=("Fujita.mpl" "Goodwin.mpl" "HIV2.mpl" "New_AKT.mpl" "OralGlucose.mpl" "QWWC.mpl" "SEIQRDC.mpl" "SEIRP.mpl" "SIRAQJ.mpl" "Treatment.mpl" "new_SEIR.mpl")
for example in "${array[@]}"; do
    /local/maple2021/bin/maple $example
done
