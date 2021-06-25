#!/bin/bash
r=100 #! require exact match up to r, within 1/r in general
s=$(grep -F ":${2}:" mzdata.txt | awk -F: '{OFS=":"; print $5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19}')
if [ "$s" == "" ]; then
    echo "ERROR: LMFDB label ${2} not found in mzdata.txt!"
    exit 99
fi
m=$(awk -F: '{OFS=":"; print $5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19}' moments_${1}_${2}.txt)
if [ "$m" == "" ]; then
    echo "ERROR: moment data not found in moments_${1}_${2}.txt"
    exit 99
fi
A=$(awk -F: '{print $2}' moments_${1}_${2}.txt)
a=(${s//:/ })
b=(${m//:/ })
if [ ${a[0]} != ${b[0]} ]; then
    echo "ERROR: moment/zvector data for Example $1 $A in moments_${2}.txt does not match $2"
    echo "Expected w2 = ${a[0]} got w2 = ${b[0]}"
    echo "Expected w4 = ${a[1]} got w4 = ${b[1]}"
    echo "Expected d10 = ${a[12]} got d10 = ${b[12]}"
    exit 99
fi
if [ ${a[1]} != ${b[1]} ]; then
    n=${#a[1]}; c=${a[1]:1:n-2}; x=(${c//,/ })
    n=${#b[1]}; c=${b[1]:1:n-2}; y=(${c//,/ })
    for i in {0..3}; do
        if [ ${x[i]} != ${y[i]} ]; then
            let "d = $r * (${x[i]} - ${y[i]})"
            if [ $d -lt 0 ]; then let "d = - d"; fi
            if [ $d -gt ${x[i]} ]; then
                echo "ERROR: moment/zvector data for Example $1 $A in moments_${2}.txt does not match $2"
                echo "Expected w4 = ${a[1]} got w4 = ${b[1]} differing in entry $((i+1))"
                echo "Expected d10 = ${a[12]} got d10 = ${b[12]}"
                exit 99
            fi
        fi
    done
fi
if [ ${a[5]} != ${b[5]} ]; then
    x=${a[5]:2:4}; y=${a[5]:2:4}
    let "d = 1$x - 1$y" 
    if [ $d -lt 0 ]; then let "d = - d"; fi
    if [ $d -gt 1 ]; then
        echo "ERROR: moment/zvector data for Example $1 $A in moments_${2}.txt does not match $2"
        echo "Expected z1 = ${a[5]} got z1 = ${b[5]}"
        echo "Expected d10 = ${a[12]} got d10 = ${b[12]}"
        exit 99
    fi
fi
if [ ${a[6]} != ${b[6]} ]; then
    n=${#a[6]}; c=${a[6]:1:n-2}; x=(${c//,/ })
    n=${#b[6]}; c=${b[6]:1:n-2}; y=(${c//,/ })
    for i in {0..4}; do
        let "d = 1${x[i]:2:4} - 1${y[i]:2:4}"
        if [ $d -lt 0 ]; then let "d = - d"; fi
        if [ $d -gt 1 ]; then
            echo "ERROR: moment/zvector data for Example $1 $A in moments_${2}.txt does not match $2"
            echo "Expected z2 = ${a[6]} got z2 = ${b[6]}"
            echo "Expected d10 = ${a[12]} got d10 = ${b[12]}"
            exit 99
        fi
    done
fi
if [ ${a[7]} != ${b[7]} ]; then
    x=${a[7]:2:4}; y=${a[7]:2:4}
    let "d = 1$x - 1$y" 
    if [ $d -lt 0 ]; then let "d = - d"; fi
    if [ $d -gt 1 ]; then
        echo "ERROR: moment/zvector data for Example $1 $A in moments_${2}.txt does not match $2"
        echo "Expected z3 = ${a[7]} got z3 = ${b[7]}"
        echo "Expected d10 = ${a[12]} got d10 = ${b[12]}"
        exit 99
    fi
fi
if [ ${a[8]} != ${b[8]} ]; then
    n=${#a[8]}; c=${a[8]:1:n-2}; x=(${c//,/ })
    n=${#b[8]}; c=${b[8]:1:n-2}; y=(${c//,/ })
    for i in {0..4}; do
        let "d = 1${x[i]:2:4} - 1${y[i]:2:4}"
        if [ $d -lt 0 ]; then let "d = - d"; fi
        if [ $d -gt 1 ]; then
            echo "ERROR: moment/zvector data for Example $1 $A in moments_${2}.txt does not match $2"
            echo "Expected z12 = ${a[8]} got z12 = ${b[8]}"
            echo "Expected d10 = ${a[12]} got d10 = ${b[12]}"
            exit 99
        fi
    done
fi
if [ ${a[9]} != ${b[9]} ]; then
    x=${a[9]:2:4}; y=${a[9]:2:4}
    let "d = 1$x - 1$y" 
    if [ $d -lt 0 ]; then let "d = - d"; fi
    if [ $d -gt 1 ]; then
        echo "ERROR: moment/zvector data for Example $1 $A in moments_${2}.txt does not match $2"
        echo "Expected z13 = ${a[9]} got z13 = ${b[9]}"
        echo "Expected d10 = ${a[12]} got d10 = ${b[12]}"
        exit 99
    fi
fi
if [ ${a[10]} != ${b[10]} ]; then
    n=${#a[10]}; c=${a[10]:1:n-2}; x=(${c//,/ })
    n=${#b[10]}; c=${b[10]:1:n-2}; y=(${c//,/ })
    for i in {0..4}; do
        let "d = 1${x[i]:2:4} - 1${y[i]:2:4}"
        if [ $d -lt 0 ]; then let "d = - d"; fi
        if [ $d -gt 1 ]; then
            echo "ERROR: moment/zvector data for Example $1 $A in moments_${2}.txt does not match $2"
            echo "Expected z23 = ${a[10]} got z23 = ${b[10]}"
            echo "Expected d10 = ${a[12]} got d10 = ${b[12]}"
            exit 99
        fi
    done
fi
if [ ${a[11]} != ${b[11]} ]; then
    n=${#a[11]}; c=${a[11]:1:n-2}; x=(${c//,/ })
    n=${#b[11]}; c=${b[11]:1:n-2}; y=(${c//,/ }) 
    for i in {0..4}; do
        let "d = 1${x[i]:2:4} - 1${y[i]:2:4}"
        if [ $d -lt 0 ]; then let "d = - d"; fi
        if [ $d -gt 1 ]; then
            echo "ERROR: moment/zvector data for Example $1 $A in moments_${2}.txt does not match $2"
            echo "Expected z123 = ${a[11]} got z123 = ${b[11]}"
            echo "Expected d10 = ${a[12]} got d10 = ${b[12]}"
            exit 99
        fi
    done
fi
# if [ ${a[14]} != ${b[14]} ]; then
#     n=${#a[14]}; c=${a[14]:1:n-2}; x=(${c//,/ })
#     n=${#b[14]}; c=${b[14]:1:n-2}; y=(${c//,/ })
#     for i in {0..6}; do
#         let "d = 1${x[i]:2:4} - 1${y[i]:2:4}"
#         if [ $d -lt 0 ]; then let "d = - d"; fi
#         if [ $d -gt 1 ]; then
#             echo "ERROR: moment/zvector data for Example $1 $A in moments_${2}.txt does not match $2"
#             echo "Expected d-vector = ${a[14]} got d-vector = ${b[14]}"
#             echo "Expected d10 = ${a[12]} got d10 = ${b[12]}"
#             exit 99
#         fi
#     done
# fi
n=${#a[2]}; c=${a[2]:1:n-2}; x=(${c//,/ })
n=${#b[2]}; c=${b[2]:1:n-2}; y=(${c//,/ })
for i in {0..6}; do
    if [ ${x[i]} != ${y[i]} ]; then
        let "d = $r * (${x[i]} - ${y[i]})"
        if [ $d -lt 0 ]; then let "d = - d"; fi
        if [ $d -gt ${x[i]} ]; then
            echo "ERROR: moment/zvector data for Example $1 $A in moments_${2}.txt does not match $2"
            echo "Expected w6 = ${a[2]} got w6 = ${b[2]} differing in entry $((i+1))"
            echo "Expected d10 = ${a[12]} got d10 = ${b[12]}"
            exit 99
        fi
    fi
done
n=${#a[3]}; c=${a[3]:1:n-2}; x=(${c//,/ })
n=${#b[3]}; c=${b[3]:1:n-2}; y=(${c//,/ })
for i in {0..6}; do
    if [ ${x[i]} != ${y[i]} ]; then
        let "d = $r * (${x[i]} - ${y[i]})"
        if [ $d -lt 0 ]; then let "d = - d"; fi
        if [ $d -gt ${x[i]} ]; then
            echo "ERROR: moment/zvector data for Example $1 $A in moments_${2}.txt does not match $2"
            echo "Expected w8 = ${a[3]} got w6 = ${b[3]} differing in entry $((i+1))"
            echo "Expected d10 = ${a[12]} got d10 = ${b[12]}"
            exit 99
        fi
    fi
done
n=${#a[4]}; c=${a[4]:1:n-2}; x=(${c//,/ })
n=${#b[4]}; c=${b[4]:1:n-2}; y=(${c//,/ })
for i in {0..6}; do
    if [ ${x[i]} != ${y[i]} ]; then
        let "d = $r * (${x[i]} - ${y[i]})"
        if [ $d -lt 0 ]; then let "d = - d"; fi
        if [ $d -gt ${x[i]} ]; then
            echo "ERROR: moment/zvector data for Example $1 $A in moments_${2}.txt does not match $2"
            echo "Expected w10 = ${a[4]} got w6 = ${b[4]} differing in entry $((i+1))"
            echo "Expected d10 = ${a[12]} got d10 = ${b[12]}"
            exit 99
        fi
    fi
done
n=${#a[12]}; c=${a[12]:1:n-2}; x=(${c//,/ })
n=${#b[12]}; c=${b[12]:1:n-2}; y=(${c//,/ })
for i in {0..9}; do
    if [ ${x[i]} != ${y[i]} ]; then
        let "d = $r * (${x[i]} - ${y[i]})"
        if [ $d -lt 0 ]; then let "d = - d"; fi
        if [ $d -gt ${x[i]} ]; then
            echo "ERROR: moment/zvector data for Example $1 $A in moments_${2}.txt does not match $2"
            echo "Expected d10 = ${a[12]} got d10 = ${b[12]} differing in entry $((i+1))"
            exit 99
        fi
    fi
done
n=${#a[13]}; c=${a[13]:1:n-2}; x=(${c//,/ })
n=${#b[13]}; c=${b[13]:1:n-2}; y=(${c//,/ })
for i in {0..9}; do
    if [ ${x[i]} != ${y[i]} ]; then
        let "d = $r * (${x[i]} - ${y[i]})"
        if [ $d -lt 0 ]; then let "d = - d"; fi
        if [ $d -gt ${x[i]} ]; then
            echo "ERROR: moment/zvector data for Example $1 $A in moments_${2}.txt does not match $2"
            echo "Expected d20 = ${a[13]} got d20 = ${b[13]} differing in entry $((i+1))"
            exit 99
        fi
    fi
done
echo "SUCCESS: moment/zvector data for for Example $1 $A matches $2"
exit 0
