######################################################
### DEFINITION OF SIGNALED TRANSPORTATION P SYSTEM ###
######################################################


### MULTISET ###
setClass('multiset'
         ,slots=c(
            obj="character",
            tar="numeric"))

# init multiset#
createMultiset = function(object=character(0),target=numeric(0)){
  mset = new('multiset',obj=object,tar=target)
  return(mset)
}

# print multiset

printMultiset = function(mset){
  if (length(mset@obj) > 0L){
    if (length(mset@tar) > 0L){
      print(paste0(mset@obj,'(in_',mset@tar,')'))
    } else {
      print(mset@obj)
    }
  } else {
    print('lambda')
  }
}

### LIST OF MULTISETS ###
setClass("multisets",slots=c(objects='list'))


### SYNAPSE BETWEEN MEMBRANES ###
setClass('synapse',slots=c(
  f="numeric",
  t="numeric",
  s="logical"))

### init synapse ###
createSynapse = function(from,to,status=logical(0)){
  syn = new('synapse',f=from,t=to,s=status)
  return(syn)
}

### MULTISET OPTIONAL - allow to be either multiset or NA
setClassUnion("synapseOpt",   c("synapse",   "logical"))


### LEFT SIDE OF RULE ###
setClass('ruleLeftSide',slots=c(
  memblabel='numeric',
  obj='multiset',
  catalyst='multiset',
  syn='synapseOpt'))

### init left side of rule ###
createRLS = function(label,object,synapse=NA,catalyst=createMultiset()){
  r <- new('ruleLeftSide',memblabel=label,obj=object,syn=synapse,catalyst=catalyst)  
  return(r)
}



### RIGHT SIDE OF RULE ###
setClass('ruleRightSide',slots=c(
  memblabel='numeric',
  obj='multiset'))

### init right side of rule ###
createRRS = function(label,object){
  r <- new('ruleRightSide',memblabel=label,obj=object)  
  return(r)
}



### RULE ###
setClass('rule',
         slots=c(
           left='ruleLeftSide',
           right='ruleRightSide',
           time='numeric',
           prob='numeric',
           statusFrom='logical',
           statusTo='logical'))

### init rule ###
createRule = function(lrs,rrs,t=1,p=1,statFrom=logical(0),statTo=logical(0)){
  r <- new('rule',left=lrs,right=rrs,time=t,prob=p,statusFrom=statFrom,statusTo=statTo)  
  return(r)
}

### print ###
printRule <- function(r){
  left <- ''
  syntext <- ''
  if (!is.na(r@left@syn)){
    syntext <- paste0('syn',r@left@syn@f,r@left@syn@t)
    if (length(r@left@syn@s)!=0L){
      if(r@left@syn@s){
        syntext <- paste0(syntext,'+')
      } else {
        syntext <- paste0(syntext,'-')
      }    
    }
  }
  catext <-''
  if (length(r@left@catalyst@obj) > 0L){
    catext <-paste0('|',r@left@catalyst@obj)
  }
  left <- paste0('[',r@left@obj@obj,r@left@obj@tar,'|',syntext,catext,']',r@left@memblabel)
  right <- ''
  right <- paste0('[',r@right@obj@obj,r@right@obj@tar,']',r@right@memblabel)
  print(paste0(left,' -> ',right,' ; ',r@time,', ',r@prob,',  ',r@statusFrom,', ',r@statusTo))
}


### Event ###
setClass('event',slots=c(rule='rule',time='numeric'))

### init event ###
createEvent = function(r,t){
  e <- new('event',rule=r,time=t)  
  return(e)
}


### LIST OF RULES ###
setClass('rules',slots=c(rules='list'))



### PENDING RULE ###
setClass('pendingRule',slots=c(
  catalyst='multiset',
  membraneFrom='numeric',
  right='ruleRightSide',
  time='numeric',
  syn='synapseOpt',
  statusFrom='logical',
  statusTo='logical'))

### init pending rule ###
createPendingRule = function(membFrom,rrs,t,synapse,statFrom=logical(0),statTo=logical(0),catalyst=createMultiset()){
  r <- new('pendingRule',catalyst=catalyst,membraneFrom=membFrom,right=rrs,syn=synapse,time=t,statusFrom=statFrom,statusTo=statTo)  
  return(r)
}

# create pending rule from regural rule
createPendingRuleFromRule = function(rule){
  r <- new('pendingRule',catalyst=rule@left@catalyst,membraneFrom=rule@left@memblabel,right=rule@right,syn=rule@left@syn,time=rule@time-1,statusFrom=rule@statusFrom,statusTo=rule@statusTo)  
  return(r)
}




### MEMBRANE ###
setClass('membrane',slots=c(
  label="numeric",
  objects="multisets",
  rules="rules"))

createMembrane <- function(memblabel,objects=new('multisets',objects=list()),rules){
  m <- new('membrane',label=memblabel,objects=objects,rules=rules)
  return(m)
}


### MEMBRANE STRUCTURE ###
setClass('membraneStructure',slots=c(
  membranes='list',
  synapses="list",
  pendingRules='list',
  pendingSynapses='list',
  events='list'))

createMembraneStructure <- function(membranes,synapses,pendingRules,pendingSynapses=list(),evts=list()){
  ms <- new('membraneStructure',
            membranes=membranes,
            synapses=synapses,
            pendingRules=pendingRules,
            pendingSynapses=pendingSynapses,
            events=evts)
}

