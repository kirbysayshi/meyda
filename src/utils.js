module.exports.µ = function(i, amplitudeSpect){
  var numerator = 0;
  var denominator = 0;
  
  for(var k = 0; k < amplitudeSpect.length; k++){
    numerator += Math.pow(k,i)*Math.abs(amplitudeSpect[k]);
    denominator += amplitudeSpect[k];
  }

  return numerator/denominator;
};

module.exports.isPowerOfTwo = function(num) {
  while (((num % 2) === 0) && num > 1) {
    num /= 2;
  }

  return (num == 1);
};