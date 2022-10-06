data {
 array[10] int counts;
}
parameters {
 simplex[10] theta;
}
model {
  counts ~ multinomial(theta);
}
