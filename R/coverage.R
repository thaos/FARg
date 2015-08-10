#coverage <- function(generate_data, get_theoric, get_ic)

#set_env set the parameters
#generate_data <- function() return data.frame/or an environment
#get_ic takes a environement, return theoric ic given the way data are generated
#get_theoric takes a environement, return CIs (inf,hat,sup) for the wanted value

coverage <- function(set_env, gen_data, get_theo, get_ic ){
  function(n, ...){
    ics <- matrix(ncol=4, nrow=n)
    colnames(ics) <- c("Theoric", "IC_inf", "Esti", "IC_sup") 
    env <- set_env()
    env <- gen_data(env)
    ics[, 1] <- get_theo(env)
    for(i in seq.int(n))
        ics[i,-1] <- get_ic(env,...)[1:3]
    #     str( vapply(seq.int(n),function(x, ...)get_ic(env, ...)[1:3],FUN.VALUE=numeric(3), ...))
        #     ics[, -1] <- t(vapply(seq.int(n),function(x,...)get_ic(env, ...)[1:3], FUN.VALUE=numeric(3), ...))
    cover <- mean((ics[,1] >= ics[,2] & ics[,1] <= ics[,4]))
    list("coverage"=cover, "IC"=ics)
  }
}
    





  
