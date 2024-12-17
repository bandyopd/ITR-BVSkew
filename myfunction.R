phi<-function(v,u1,u2){
  return( v*(1-u1-u2+v)/((u1-v)*(u2-v)) )
}

inv_phi<-function(v,u1,u2){
  a=v-1
  b=u1+u2-(u1+u2)*v-1
  c=v*u1*u2
  output=(-b-sqrt(pmax(b^2-4*a*c,0)   ))/(2*a)
  if(sum(v==1)>0){
    output[v==1]=(u1*u2)[v==1]
  }
  output[is.na(output)]=runif(sum(is.na(output)))
  return(output)
}


####### objective function for bivariate distributional regression
ass_f<-function(t,w){
  R_a=exp(tempX%*%t)
  prob=inv_phi(R_a,new1,new2)
  value=-sum((1/w)* ( O_L*log(pmax(prob,1e-3) )+
               (1-O_L)*log(pmax(1-prob,1e-3)) ) )
  return(value)
}
