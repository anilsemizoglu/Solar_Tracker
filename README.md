A solar tracker developed in Visual Studio 2012, using a Mightex CCD USB camera.

#Purpose and Tools

The purpose of this piece of software is to provide information to a telescope's motors about which way to aim the telescope.
Ephemeris file used has the precise coordinates of the sun in the astronomical coordinates.
If the camera and a telescope is aligned so they're both on the same beam axis, the centroid location of the sun in a captured frame corresponds to the real location of the sun in the space, hence we can adjust the telescope to look at a specific point in the space.

#Libraries and Dependencies:

  + JPL ephemeris ASCII format (using 406 in this version, can use an updated version although not tested.)- ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/

  + JPL ephemeris ASCII to binary conversion (asc2bin.c)- http://www.cv.nrao.edu/~rfisher/Ephemerides/computer_code.html#asc2bin

  + SofaJPL (to extraxt information from the JPL ephemeris) - http://home.earthlink.net/~s543t-24dst/sofajpl/
  
  + OpenCV 2.4.11 (image processing)- http://sourceforge.net/projects/opencvlibrary/files/opencv-win/2.4.11/
  
  + Link the libraries in project settings in Visual Studio.
 
#Source Files:

Solar_Tracker.cpp:
 
	+ main function is located in this file.
	
	+ Settings are loaded from the Settings.txt
 	
 	+ FrameCallBack is where the frame capture and image processing takes place.
 
 	+ outputs a CSV file of the timestamp, apparent RA/DEC position of the sun at the 
 	time of the frame capture, and the centroid location of the sun in the image frame in pixel numbers.
 
 eph_calc.cpp:
 
 	+ has the sun_equatorial function that calculates the RA/DEC position of the sun.
 	
 	+ loads the binary JPL ephemeris ( "binp2000.406" in this case )
 	
utility.cpp:
 	
 	+ get_date function for the timestamp, get_cpu_time for performance checks.
 	
ini_reader.cpp:
 	
 	+ search functions to load the settings file.
