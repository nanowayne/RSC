void findradialsymmetrycenter(image select, number &xc, number &yc)
{
// For information on this method see Raghuveer Parthasarathy - Nature Methods 9, 724–726 (2012)
// This program adpate from his original matlab code
// http://physics-server.uoregon.edu/~raghu/Particle_tracking_files/radialcenter.m

number xsize, ysize
getsize(select, xsize, ysize)

// This what i understancd 
// and nomornize the weight
Image weightx = RealImage("", 4, xsize, ysize)
weightx =-(xsize-2)/2+icol
Image weighty = RealImage("", 4, xsize, ysize)
weighty =-(ysize-2)/2+irow

// Calculate derivatives along 45-degree shifted coordinates (u and v)
//By offset commend
Image dIdu = CreateFloatImage("dIdu",xsize,ysize)
dIdu= offset(select,1,0)-offset(select,0,1)
Image dIdv = CreateFloatImage("dIdv",xsize,ysize)
dIdv= offset(select,0,0)-offset(select,1,1)

//smoothing by gaussian kernel
Image gaussiankernel:= [3,3]:
			{
				{1,1,1},
				{1,1,1},
				{1,1,1}
			}
			
Image fdu = Convolution(dIdu, gaussiankernel/9)	
Image fdv = Convolution(dIdv, gaussiankernel/9)	
Image dImg2 = fdu*fdu+fdv*fdv

//slope of the gradient
Image m = - (fdv+fdu)/(fdu-fdv)


//If the m is Nan and infinite
m = tert(IsNan(m), 0, m)
m = tert(IsInfinite(m), 0, m)

Image b = weighty - m*weightx

// Weighting: weight by squre of gradient magnitude and inverse
number xcentroid, ycentroid

xcentroid = sum(dImg2*weightx)/sum(dImg2)
ycentroid = sum(dImg2*weighty)/sum(dImg2)
Image w
w = dImg2/sqrt((weightx-xcentroid)*(weightx-xcentroid)+(weighty-ycentroid)*(weighty-ycentroid))

// Least-square minization for center
number sw, smmw, smw, smbw, sbw, det
Image wm2pl = w/(m*m +1)
sw = sum(wm2pl)
smmw = sum(m*m*wm2pl)
smw = sum(m*wm2pl)
smbw = sum(m*b*wm2pl)
sbw = sum(b*wm2pl)
det = smw*smw - smmw*sw
xc = (smbw*sw - smw*sbw)/det //relative to the image center
yc = (smbw*smw - smmw*sbw)/det //reolative to the image center

//retun the output relate to the upper left coordinates
xc = xc + (xsize+1)/2-1
yc = yc + (ysize+1)/2-1

}
  
// Main script 
// have an image front-most which peak (atom column) feature, which you
// wish to locate the atomic column center by radial symmetry center method
 
// Check the time
number ticks=GetHighResTickCount()

image temp:=getfrontimage()
number xc, yc
findradialsymmetrycenter(temp, xc,yc)
result("\nCx="+xc+" Cy="+yc)

Createovalannotation(temp, yc-5, xc-5, yc+5, xc+5) 

// Display the elapsed time
number lastticks=GetHighResTickCount()
number elapsed=CalcHighResSecondsBetween(ticks, lastticks)
result("\n\nElapsed Time = "+elapsed+" s\n")