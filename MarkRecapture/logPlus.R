log.plus<- function(x,y) 
{
  # Given x=log(a) and y=log(b), return log(a+b)
  # Trick by John Skilling.
  if(x>y) 
  {
    x+log(1+exp(y-x))
  } else
  {
    y+log(1+exp(x-y))
  }
}
