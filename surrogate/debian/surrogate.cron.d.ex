#
# Regular cron jobs for the surrogate package
#
0 4	* * *	root	[ -x /usr/bin/surrogate_maintenance ] && /usr/bin/surrogate_maintenance
