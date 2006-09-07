require("fracdiff")

set.seed(17)
memory.long <- fracdiff.sim(1500, d = 0.3)
spm <- fdSperio(memory.long$series)
str(spm, digits=6)

whittle1(memory.long$series)
whittle2(memory.long$series)

## Plot the *local* whittle  "likelihood: -- but not just in  [0.5, 1] :
curve(sapply(x, function(d) lw(d, memory.long$series, 5)), -.1, 1)
