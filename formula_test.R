data(tas)
f <- eur_tas ~ avg_gbl_tas + bob
f <- eur_tas ~ 1
f <- ~ avg_gbl_tas
f <- ~ 
class(f)
terms(f)
model.matrix(f, data=tas)
complete_formula("x", f, env=environment()) 
x <- 1
complete_formula( x, f, env=environment()) 
complete_formula( y, f)
lm(complete_formula( "eur_tas", f),data=tas)
complete_formula( "eur_tas", f)

