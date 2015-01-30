module.exports = function(bufferSize, m) {
  var ampspec = m.ampSpectrum;
  var µ1 = µ(1, ampspec);
  var µ2 = µ(2, ampspec);
  var µ3 = µ(3, ampspec);
  var µ4 = µ(4, ampspec);
  var numerator = -3 * Math.pow(µ1, 4) + 6 * µ1 * µ2 - 4 * µ1 * µ3 + µ4;
  var denominator = Math.pow(Math.sqrt(µ2 - Math.pow(µ1, 2)), 4);
  return numerator / denominator;
};