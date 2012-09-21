# Set the task name
TASK = fss_check

# Uncomment the correct choice indicating either SKA or TST flight environment
FLIGHT_ENV = SKA

# Set the names of all files that get installed
SHARE = daily_fss.py check_fss.py bad_times.py
DATA = task_schedule.cfg
WWW = index.html

include /proj/sot/ska/include/Makefile.FLIGHT
INSTALL_WWW = /data/mta4/www/ASPECT/$(TASK)

install:
#  Uncomment the lines which apply for this task
	mkdir -p $(INSTALL_SHARE)
	mkdir -p $(INSTALL_DATA)/fssa
	mkdir -p $(INSTALL_DATA)/fssb
	mkdir -p $(INSTALL_WWW)/fssa
	mkdir -p $(INSTALL_WWW)/fssb
	rsync --times --cvs-exclude $(SHARE) $(INSTALL_SHARE)/
	rsync --times --cvs-exclude $(DATA) $(INSTALL_DATA)/
	rsync --times $(WWW) $(INSTALL_WWW)
