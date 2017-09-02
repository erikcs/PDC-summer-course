library(ggplot2)

# First parallel version
#t = c(109, 73, 77, 103, 96, 81)
t = c(60, 70, 77, 102, 96, 81)
t = t/60
c = c(1, 2, 4, 8, 16, 32)
nobs = length(t)
v = replicate(nobs, "1st (inner loop, race condition)")
g = replicate(nobs, 5000)
s = replicate(nobs, "default")
i = replicate(nobs, 1)
x = seq(1, nobs)
df1 = data.frame(
  v, c, t, g, s, i, x
)

#Second parallel (outer loop), n = 5000 Note: ran these again to have comparable timings for the same day...
# old: t = c(54, 33, 17, 9, 5, 2)
t = c(47, 26.7, 14.2, 8, 4.5, 2.5 )
t = t/47
c = c(1, 2, 4, 8, 16, 32)
nobs = length(t)
v = replicate(nobs, "v2 outer loop, static schedule") ## "Outer loop pragma") #"v2 (outer loop, static schedule)")
g = replicate(nobs, 5000)
s = replicate(nobs, "default")
i = replicate(nobs, 2)
x = seq(1, nobs)
df2 = data.frame(
  v, c, t, g, s, i, x
)

#Second parallel (outer loop), n = 8000
t = c(142, 94, 44, 23, 12, 6) 
t = t/142
c = c(1, 2, 4, 8, 16, 32)
nobs = length(t)
v = replicate(nobs, "Second parallel (preceding loop)")
g = replicate(nobs, 8000)
s = replicate(nobs, "default")
i = replicate(nobs, 3)
x = seq(1, nobs)
df3 = data.frame(
  v, c, t, g, s, i, x
)

#Second parallel (outer loop), n = 100
t = c(0.007, 0.08)
c = c(1, 32)
nobs = length(t)
v = replicate(nobs, "Second parallel (preceding loop)")
g = replicate(nobs, 100)
s = replicate(nobs, "default")
i = replicate(nobs, 4)
x = seq(1, nobs)
df4 = data.frame(
  v, c, t, g, s, i, x
)

#Second parallel (outer loop), n = 5000, dynamic scheduling
t = c(46.8, 78.3, 85.3, 75.1, 57.3, 40.4)
t = t/46.8
c = c(1, 2, 4, 8, 16, 32)
nobs = length(t)
v = replicate(nobs, "v2 (outer loop, dynamic scheduling)")
g = replicate(nobs, 5000)
s = replicate(nobs, "dynamic")
i = replicate(nobs, 5)
x = seq(1, nobs)
df5 = data.frame(
  v, c, t, g, s, i, x
)

#Second parallel (outer loop), n = 5000, guided scheduling
t = c(46.8, 24, 13, 8, 6.2, 5.6)
t = t/46.8
c = c(1, 2, 4, 8, 16, 32)
nobs = length(t)
v = replicate(nobs, "v2 (outer loop, guided scheduling)")
g = replicate(nobs, 5000)
s = replicate(nobs, "guided")
i = replicate(nobs, 6)
x = seq(1, nobs)
df6 = data.frame(
  v, c, t, g, s, i, x
)

#Th. runtime....
t = c(1, 1/2, 1/2^2, 1/2^3, 1/2^4, 1/2^5)
c = c(1, 2, 4, 8, 16, 32)
nobs = length(t)
v = replicate(nobs, "Theoretical (1/2^threads)")
g = replicate(nobs, 5000)
s = replicate(nobs, "guided")
i = replicate(nobs, 7)
x = seq(1, nobs)
df7 = data.frame(
  v, c, t, g, s, i, x
)

df = rbind(df1, df2, df3, df4, df5, df6, df7)
names(df) = c("version", "threads", "runtime", "grid", "tscheduling", "id", "x")

p = ggplot(df[df$grid==5000, ],
           aes(x=x, runtime, color=factor(version))) + 
             geom_point() + geom_line()  +
              scale_x_discrete(labels=c(1, 2, 4, 8, 16, 32)) +  
                labs(x="threads", y="runtime (sec)") + theme(legend.title=element_blank()) +
                  ggtitle("Parallel runtimes for grid size = 5000") +
  scale_y_continuous(breaks = round(seq(min(df$runtime), max(df$runtime), by = .5),1))
 
p
setwd("/Users/erik/Dropbox/Akademisk/DN2258 PDC/Project/latex")
#ggsave("plot1.pdf")

#Last: combined 2 loops, parallel chunck
t = c(94, 55, 29, 15, 8, 4)
t = t/94
c = c(1, 2, 4, 8, 16, 32)
nobs = length(t)
v = replicate(nobs, "Combined loops")
g = replicate(nobs, 5120)
s = replicate(nobs, "static")
i = replicate(nobs, 6)
x = seq(1, nobs)
dfLast = data.frame(
  v, c, t, g, s, i, x
)

#Last: inner...
t = c(56, 38, 27,24 , 24, 26)
t = t/56
c = c(1, 2, 4, 8, 16, 32)
nobs = length(t)
v = replicate(nobs, "Inner loop pragma")
g = replicate(nobs, 5120)
s = replicate(nobs, "static")
i = replicate(nobs, 6)
x = seq(1, nobs)
dfLast2 = data.frame(
  v, c, t, g, s, i, x
)

dfLast = rbind(dfLast, dfLast2, df2)
names(dfLast) = c("version", "threads", "runtime", "grid", "tscheduling", "id", "x")


####
p01 = ggplot(df[df$i==1 , ],
           aes(x=x, runtime) )+ 
  geom_point() + geom_line()  +
  scale_x_discrete(labels=c(1, 2, 4, 8, 16, 32)) +  
  labs(x="threads", y="runtime") + theme(legend.title=element_blank()) +
  ggtitle("Runtime (1 thread = 100%), Grid size = 5000, first version (inner loop parallization)") +
  scale_y_continuous(breaks = round(seq(min(df$runtime), max(df$runtime), by = .2),1))

p01

p02 = ggplot(df[df$grid==5000 & df$i!=1, ],
           aes(x=x, runtime, color=factor(version))) + 
  geom_point() + geom_line()  +
  scale_x_discrete(labels=c(1, 2, 4, 8, 16, 32)) +  
  labs(x="threads", y="runtime") + theme(legend.title=element_blank()) +
  ggtitle("Runtime (1 thread = 100%), Grid size = 5000, second version (outer loop)") +
  scale_y_continuous(breaks = round(seq(min(df$runtime), max(df$runtime), by = .5),1))

p02

p03 = ggplot(dfLast,
            aes(x=x, runtime, color=factor(version))) + 
  geom_point() + geom_line()  +
  scale_x_discrete(labels=c(1, 2, 4, 8, 16, 32)) +  
  labs(x="threads", y="runtime") + theme(legend.title=element_blank()) +
   
  scale_y_continuous(breaks = round(seq(min(df$runtime), max(df$runtime), by = .1),1))

p03

#Last: inner...
t = c(3.69, 4.13, 4.36,4.56 , 4.89, 5.07)
t = t/3.69
c = c(1, 2, 4, 8, 16, 32)
nobs = length(t)
v = replicate(nobs, "Inner loop pragma")
g = replicate(nobs, 5120)
s = replicate(nobs, "static")
i = replicate(nobs, 6)
x = seq(1, nobs)
dfLast3 = data.frame(
  v, c, t, g, s, i, x
)
names(dfLast3) = c("version", "threads", "runtime", "grid", "tscheduling", "id", "x")
p04 = ggplot(dfLast3,
             aes(x=x, runtime)) + 
  geom_point() + geom_line()  +
  scale_x_discrete(labels=c(1, 2, 4, 8, 16, 32)) +  
  labs(x="threads", y="runtime") + theme(legend.title=element_blank()) +
  ggtitle("Weak scaling") +
  scale_y_continuous(breaks = round(seq(0.5, 1.5, by = .1),1)) + ylim(0.5, 1.5)

p04




