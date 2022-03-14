array=("Goodwin.mpl" "OralGlucose.mpl" "QWWC.mpl" "SEIAJRC.mpl" "SIAR.mpl" "SEIRP.mpl" "SIRAQJ.mpl")
for example in "${array[@]}"; do
    maple $example
done
