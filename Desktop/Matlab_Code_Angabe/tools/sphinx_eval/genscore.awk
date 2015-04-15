#=======================================================
#
# Produce keyword (letters and digits) recognition accuracy
# from CHiME result files in the following format:
#
#     s10_bgakzn k z
#     s10_bgay1a y 9
#     s10_bgiq3a u 3
#
# Ning Ma, Sheffield University
# 12 Jan 2010
#=======================================================

BEGIN {
  cnt = 0; score = 0;
}
{
  # Letter
  if (substr($1,4,1)==$2) score++;
  # Digit
  if (substr($1,5,1)==$3) score++; 

  cnt++;
}
END { 
  # Only output results if all the utterances are finished
  printf("%.2f\n", score*100/2/cnt);
}

