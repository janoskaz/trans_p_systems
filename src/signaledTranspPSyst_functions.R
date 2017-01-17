
########################################
### FUNCTIONS TO CHECK APPLICABILITY ###
########################################

# check, if membrane exists - deprecated
# IN: membraneStructure, numeric (membrane label)
# OUT: logical
# membrane structure must be declared
checkMembraneExists <- function(label){
  return(memb %in% unlist(
    lapply(
      ms@membranes,
      FUN=function(m){return(m@label)}
    )))
}

# check, if synapse s2 is equal to synapse s1
# IN: synapse, synapse
# OUT: logical
compareSynapses <- function(s1,s2){
  # labels are different
  if (s1@f != s2@f | s1@t != s2@t){return(FALSE)}
  # if status is set, compare them
  if (length(s1@s) != 0L){
    return (s1@s == s2@s)
  } else {
    return(TRUE)
  }
}

# check, if synapse exists
# IN: membraneStructure, synapse
# OUT: logical
# membrane structure must be declared
checkSynapseExists <- function(syn,synapselist){
  return(
    any(
      unlist(
        lapply(
          synapselist,
          FUN=function(x){
            return(compareSynapses(syn,x))
          }
        )
      )))
}

# check, if synapse s2 is equal to synapse s1
# IN: multiset, multiset
# OUT: logical
compareObjects <- function(o1,o2){
  # one of objects is empty
  if (length(o1@obj) == 0L | length(o2@obj) == 0L){return(FALSE)}
  # labels are different
  if (o1@obj != o2@obj){return(FALSE)}
  # if status is set, compare them
  if(length(o1@tar)!=0L){
    if(length(o2@tar)!=0L){
      return (o1@tar == o2@tar)
    } else {
      return(FALSE) 
    }    
  } else {
    # are both targets not set?
    return(length(o2@tar)==0L)
  }
}

# check, if synapse exists
# IN: multiset, synapse
# OUT: logical
checkObjectsPresent <- function(mset,obj){
  return(
    any(unlist(
      lapply(
        mset,
        FUN=function(x){
          # object is not set or objects are equal
          return(length(obj@obj) == 0L | compareObjects(obj,x))
        }
      )
    )))
}

# check, if synapse can be changed
checkSynapseCanBeChanged <- function(rule){
  if (length(rule@statusFrom) == 0L | length(rule@statusTo) == 0L){
    return(TRUE)
  } else {
    newsyn <- createSynapse(rule@left@syn@f,rule@left@syn@t)
    return(!checkSynapseExists(newsyn,ms@pendingSynapses))
  }
}

# is rule applicable
# IN: rule, membrane (current membrane, in which the rule is checked)
# OUT: logical
# membrane structure must be declared as ms
checkApplicability <- function(rule,currentMembrane){
  # objects
  left <- rule@left
  catalyst <- rule@left@catalyst
  memb <- left@memblabel
  obj <- left@obj
  syn <- left@syn
  # membrane exists
  #if (!checkMembraneExists(ms,memb)){return(FALSE)}
  # synapse exists
  if(!is.na(syn)){
    if (!checkSynapseExists(syn,ms@synapses)){return(FALSE)}
  }
  # if synapse should be changed, check pending rules
  if(!checkSynapseCanBeChanged(rule)){return(FALSE)}
  # if catalyst is present, check presence in membrane
  if (length(catalyst)!=0){
    if(!checkObjectsPresent(currentMembrane@objects@objects,catalyst)){
      return(FALSE)
    }
  }
  # objects present
  if (length(obj@obj)!=0){
    # non-empty multiset
    return( checkObjectsPresent(currentMembrane@objects@objects,obj) ) 
  } else {
    # object not present
    return(TRUE)
  }
}

###########################################
### SELECT APPLICABLE RULES IN MEMBRANE ###
###########################################

# selects all aplicable rules in membrane
# IN: membrane (curent membrane)
# OUT: list (rules)
# membrane structure must be declared
selectApplicableRules <- function(currentMembrane){
  
  lst <- lapply(
    currentMembrane@rules@rules,
    FUN=function(x){
      if(checkApplicability(x,currentMembrane)){return(x)}
      })
  # remove null
  filtered <- Filter(Negate(function(x) is.null(unlist(x))), lst)
  return(filtered)
}

# chooses rule randomly based on probabilities
# IN: list (rules)
# OUT: numeric (rule position)
# membrane structure must be declared as ms
chooseRule <- function(rules){
  # get (pseudo)probabilities
  probs <- unlist(lapply(rules,function(x){x@prob}))
  # ramdomly choose rule position
  rulePosition <- sample(length(rules),size=1,prob=probs)
  return(rulePosition)
}

##################
### APPLY RULE ###
##################

# insert pending synapse into list of pending synapses
# IN: rule
# OUT: void
# membrane structure must be declared as ms
insertPendingSynapse <- function(rule){
  if (length(rule@statusFrom)!=0L){
    newsyn <- rule@left@syn
    newsyn@s = logical(0)
    ms@pendingSynapses[[length(ms@pendingSynapses)+1]] <<- newsyn
  }
}

# apply rule
# IN: rule
# OUT: void
# membrane structure must be declared as ms
applyRule <- function(rule){
  rulesApplied <<- TRUE
  # find position of m,embrane
  label <- rule@left@memblabel
  label2 <- rule@right@memblabel
  wm <- which(unlist(lapply(ms@membranes,FUN=function(x)return(x@label == label))))
  # find position of object in membrane - only if some objects are present in membrane
  if (length(ms@membranes[[wm]]@objects@objects) > 0L){
    wo <- which(
      unlist(
        lapply(
          ms@membranes[[wm]]@objects@objects,
          FUN=function(x)return(compareObjects(rule@left@obj,x)))
      ))
    # in case there are more of the same objects
    if (length(wo)>1){
      wo <- wo[1]
    }
    # remove object
    if (length(wo) > 0){
      ms@membranes[[wm]]@objects@objects[[wo]] <<- NULL
    }
  }
  
  # remove catalysts, if necessary
  if (length(rule@left@catalyst@obj)!=0){
    wc <- which(
      unlist(
        lapply(
          ms@membranes[[wm]]@objects@objects,
          FUN=function(x)return(compareObjects(rule@left@catalyst,x)))
      ))
    # in case there are more of the same objects
    if (length(wc)>1){
      wc <- wc[1]
    }
    # remove object
    ms@membranes[[wm]]@objects@objects[[wc]] <<- NULL
  }
  # time and pending rules
  #if (rule@time>1){
  # create pending rule
  newrule <- createPendingRuleFromRule(rule)
  ms@pendingRules[[length(ms@pendingRules)+1]] <<- newrule
  # insert pending synapse, if necessary
  insertPendingSynapse(rule)
  #} else {
  #  # introduce object from right side to respective membrane
  #  wm2 <- which(unlist(lapply(ms@membranes,FUN=function(x)return(x@label == label2))))
  #  ms@membranes[[wm2]]@objects@objects <<- c(ms@membranes[[wm2]]@objects@objects, rule@right@obj)
  #}
}

# apply pending rule
# IN: pendingRule
# OUT: rule
# membrane structure must be declared as ms
applyPendingRule <- function(rule){
  rulesApplied <<- TRUE
  if (length(rule@catalyst@obj)!=0){
    wk <- which(unlist(lapply(ms@membranes,FUN=function(x)return(x@label == rule@membraneFrom))))
    ms@membranes[[wk]]@objects@objects <<- c(ms@membranes[[wk]]@objects@objects, rule@catalyst)
    rule@catalyst = createMultiset()
  }
  if (rule@time > 0){
    rule@time <- rule@time - 1
    return(rule)
  } else {
    # introduce objects into target membrane
    label <- rule@right@memblabel
    if (length(rule@right@obj@obj) > 0L){
      wm <- which(unlist(lapply(ms@membranes,FUN=function(x)return(x@label == label))))
      ms@membranes[[wm]]@objects@objects <<- c(ms@membranes[[wm]]@objects@objects, rule@right@obj)
    }    
    # change synapse state, if possible
    if (length(rule@statusFrom) !=0L & length(rule@statusTo) !=0L){
      #print(paste('changing synapse state',rule@syn@f, rule@syn@t,'; new status is',rule@statusTo))
      # change synapse state
      ws <- which(unlist(lapply(ms@synapses,FUN=function(x)return(x@f == rule@syn@f & x@t == rule@syn@t))))
      ms@synapses[[ws]]@s <<- rule@statusTo
      # remove pending synapse
      if(length(ms@pendingSynapses) != 0L){
        wps <- which(unlist(lapply(ms@pendingSynapses,FUN=function(x){return(x@f == rule@syn@f & x@t == rule@syn@t)})))
        ms@pendingSynapses[[wps[1]]] <<- NULL
      }      
    }
    # return false to show
    return(FALSE)
  }
}

# apply pending rules - all in membrane system
# IN: void
# OUT: void
# membrane structure must be declared as ms
applyPendingRules <- function(){
  rmPos = c()
  if(length(ms@pendingRules)>0){
    for (i in 1:length(ms@pendingRules)){
      res <- applyPendingRule(ms@pendingRules[[i]])
      if (class(res) == 'pendingRule'){
        ms@pendingRules[[i]] <<- res
      } else {
        rmPos <- c(rmPos,i)
      }
    }
    if (length(rmPos)>0L){
      for (j in sort(rmPos,decreasing=T)){
        ms@pendingRules[[j]] <<- NULL
      }
    }
  }  
}

# check and apply events
# IN: numeric
# OUT: void
# membrane structure must be declared as ms, indicator of rule application as rulesApplied
applyEvents <- function(t){
  for (evt in ms@events){
    if (evt@time == t){
      applyRule(evt@rule)
    }
  }
}


########################################
### FUNCTIONS TO PERFORM COMPUTATION ###
########################################

# perform one computational step in membrane
# IN: membrane
# OUT: void
# membrane structure must be declared as ms
computeInMembrane <- function(membraneIndex){
  # all applicable rules
  applicableRules <- selectApplicableRules(ms@membranes[[membraneIndex]])
  # repeat while rules can be applied
  while(length(applicableRules)>0){
    # select rule
    selectedRuleIndex <- chooseRule(applicableRules)
    #Rule(applicableRules[[selectedRuleIndex]])
    # apply rule
    applyRule(rule=applicableRules[[selectedRuleIndex]])
    # re-calculate list of applicable rules
    applicableRules <- selectApplicableRules(ms@membranes[[membraneIndex]])
  }
}

# computational step in membrane
# IN: void
# OUT: void
# membrane structure must be declared as ms
computationalStep <- function(){
  for (i in 1:length(ms@membranes)){
    computeInMembrane(i)
  }
  applyPendingRules()
}

# compute the whole system
# IN: membraneSystem, numeric
# OUT: list (of configurations in time steps)
# membrane structure must be declared as ms
simulate <- function(membraneSystem,time=Inf,fullConfig=TRUE){
  ms <<- membraneSystem
  output <- list()
  ti <- 1
  print('--- starting computation ---')
  while (ti <= time){
    rulesApplied <<- FALSE
    print(paste('step',ti,'of',time))
    applyEvents(ti)
    computationalStep()
    if (fullConfig){
      output[[ti]] <- getFullConfiguration(ms)
    } else {
      output[[ti]] <- getConfiguration(ms)
    }  
    # TODO - save configuration
    if(!rulesApplied){
      print('no rules can be applied, exit')
      return(output)
    }
    ti <- ti+1    
  }
  return(output)
}


################################
### GET SYSTEM CONFIGURATION ###
################################

# is object vehicle? regardless to target
# IN: multiset
# OUT: logical
# membrane structure must be declared as ms
isVehicle <- function(obj){
  if (obj@obj == 'veh'){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# number of vehicles in membrane
# IN: membrane
# OUT: numeric
# membrane structure must be declared as ms
getVehicleNr <- function(membrane){
  return(
    sum(
      unlist(
        lapply(
          membrane@objects@objects,
          function(x){return(isVehicle(x))}
          ))))
}

# number of vehicles in membrane
# IN: membraneSystem
# OUT: list
# membrane structure must be declared as ms
getConfiguration <- function(membraneSystem){
  return( 
    matrix(
      unlist(
        lapply(
          membraneSystem@membranes,
          function(x){return(c(x@label,getVehicleNr(x)))}
    )),nrow=2))
}

# all objects in membrane
# IN: membraneSystem
# OUT: list
# membrane structure must be declared as ms
getFullConfiguration <- function(membraneSystem){
  return( 
        lapply(
          membraneSystem@membranes,
          function(x){return(list(x@label,x@objects@objects))})
  )
}