/*
    Define your potential here as:
    
    double V(double x){
        ... do stuff
        return the value of the potential at position x.
    }
*/

//----------------------------------------------------------------
//----------------------- POTENTIAL ------------------------------
//----------------------------------------------------------------
/*
  Here you define the potential of the system V(x)
 */
double V(double x){
  double lambda = 0.1;
  return x*x + lambda*x*x*x*x;
}