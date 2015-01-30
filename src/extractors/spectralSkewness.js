module.exports = function(bufferSize, m, spectrum) {
  var ampspec = m.ampSpectrum;
  var µ1 = µ(1, ampspec);
  var µ2 = µ(2, ampspec);
  var µ3 = µ(3, ampspec);
  var numerator = 2 * Math.pow(µ1, 3) - 3 * µ1 * µ2 + µ3;
  var denominator = Math.pow(Math.sqrt(µ2 - Math.pow(µ1, 2)), 3);
  return numerator / denominator;
};