#install.packages("JuliaCall") #Do this once...

library(JuliaCall)
#julia_setup(installJulia = TRUE)

#Do this once after installing Julia 1.5
#julia_setup(JULIA_HOME = "/Applications/Julia-1.5.app/Contents/Resources/julia/bin")

getwd()
setwd("~/git/mine/2021SSA-Julia/")

julia_command( "println(\"Hello World!\")"  )


julia_command("include(\"Minimal-Simulation.jl\");")

julia_command("init_params();")
for(i in 1:10){
  results <- julia_eval("run_sim()",need_return="R")
  plot(results$time_range,results$counts_I)
}
