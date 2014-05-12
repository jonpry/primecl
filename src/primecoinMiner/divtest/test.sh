diff <(perl -pe "s/05: 0*//; tr/A-Z/a-z/" <(./a.out | grep 05:)) <(perl -pe "s/I: //" <(./a.out | grep I:))

