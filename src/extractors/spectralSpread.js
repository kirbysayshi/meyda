module.exports = function(bufferSize, m) {
  var ampspec = m.ampSpectrum;
  return Math.sqrt(µ(2, ampspec) - Math.pow(µ(1, ampspec), 2));
};