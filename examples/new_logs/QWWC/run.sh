rm *.out && touch time.out
for each in $(ls *.mpl); do
    echo $each >>time.out
    (/usr/bin/time -f "elapsed,user+system,memory\n%e,=%U+%S,%M000\n" /local/maple2021/bin/maple $each) >>output.out 2>>time.out
done
