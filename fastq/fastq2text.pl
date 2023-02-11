#print only sequences which is second line in each entry
	# @FS10000125:22:BPA73
	# GGGTTGGGCCTTTGTTTCTT
	# +
	# ,:FFFFFFF,FFFF,,,F:,

awk 'NR%4 == 2' # this will do it

# or Perl one liner
perl -ne 'print if ($.%4==2)'

