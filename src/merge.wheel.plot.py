import cv2
import 	argparse
import	cv2
import 	numpy as np
from	PIL import Image, ImageSequence
import	os
import 	re
import 	itertools
from 	matplotlib import pyplot as plt 
from skimage import data, exposure, img_as_float
import numpy as np
import sys
img=cv2.imread(sys.argv[1]+'/wheel.jpg',-1)
h,w,c=img.shape
img2=img.copy()
#img2[780:783,540:543,]=0
#cv2.imwrite('0.tif',img2[:,0:int(w*4/5),:])
#wl,wr,ht,he=[min(np.where(wheel<255)[1]),max(np.where(wheel<255)[1]),min(np.where(wheel<255)[0]),max(np.where(wheel<255)[0])]
region=img[153:1034,137:1018,]
#cv2.imwrite('1.jpg',img[153:1034,137:1018,])
### wheel[93:505,68:486,]   wheel region
font= cv2.imread('/export/home/kongyan/software/cancer_hall_mark_enrichment/database/font.jpg')
##cv2.imwrite('1.jpg',font[665:2000,1005:2340,])
## font region
## font[665:2000,1005:2340,]
fold_large=np.float32(1316/881)
enlarge = cv2.resize(region, (0, 0), fx=fold_large, fy=fold_large, interpolation=cv2.INTER_CUBIC)  
#(_, thresh) = cv2.threshold((255-enlarge),0,1, cv2.THRESH_BINARY)
#val=np.multiply(thresh,np.uint8(enlarge))
#cv2.imwrite('0.tif',enlarge)
newimg=font.copy()
newimg[665:1981,1007:2323,]=enlarge+newimg[665:1981,1007:2323,]-255
#cv2.imwrite('1.jpg',newimg)

#cv2.imwrite('1.jpg',newimg[:,3000:35000,])

#cv2.imwrite('1.jpg',img[500:700,1150:1200,])
legend=img[550:700,1150:1200,]
enlarge = cv2.resize(legend, (0, 0), fx=4, fy=4, interpolation=cv2.INTER_CUBIC) 

newimg[1750:2350,3000:3200]=enlarge
cv2.imwrite(sys.argv[1]+'/cancer_hallmark_enrichment.jpg',newimg[0:2500,0:3500,])
def save_img(img,path,name):
	cv2.imwrite(path+name,img)
 

def histeq(im,nbr_bins):
	imhist,bins = np.histogram(im.flatten(),nbr_bins,normed=True)
	# crop into nbr_bins with hists of every bins
	cdf = imhist.cumsum() # cumulative distribution function
	cdf = nbr_bins * cdf / cdf[-1]  
	#bins=nbr_bins *bins/bins[-1]
	im2 = np.interp(im.flatten(),bins[:-1],cdf)
	return im2.reshape(im.shape)

def histeq3(im,nbr_bins):
	imhist,bins = np.histogram(im.flatten(),nbr_bins,normed=False)
	cdf = imhist.cumsum() # cumulative distribution function
	cdf = nbr_bins * cdf / cdf[-1]  
	im2 = np.interp(im.flatten(),bins[:-1],cdf)
	return im2.reshape(im.shape)


def histeq2(im,nbr_bins):
	im2=np.float32(im-im.min())*np.float32(nbr_bins)/np.float32(im.max()-im.min())
	return im2

	
def cross_start_end(arr,pos,end):
	l=len(arr)
	pos2=[]
	for i in range(l):
		if arr[i]<=end:
			pos2.append(pos[i])
	return(pos2)

def mark_pos(big_figure,kernel,gray):
	match = cv2.matchTemplate(big_figure, kernel, cv2.TM_SQDIFF)
	pos=[]
	pos1=np.where(match==0)[0]
	pos2=np.where(match==0)[1]
	pos3=pos1[pos1>20] 
	pos4=pos2[pos1>20] 
	pos1=pos3[pos4>20] 
	pos2=pos4[pos4>20]
	pos3=pos1[pos1<gray.shape[0]-20]
	pos4=pos2[pos1<gray.shape[0]-20]
	pos1=pos3[pos4<gray.shape[1]-20]
	pos2=pos4[pos4<gray.shape[1]-20]
	pos=[]
	for i in range(len(pos1)):
		if max_near(gray,[pos1[i]+2,pos2[i]+4])[0]>100:
			pos.append((pos1[i]+2,pos2[i]+4))
	 
	return(pos)

def pixel_scale(img,start,end,scale_arg):
	img2=img.copy()
	img2[cross_start_end(img2[img2>=start],np.where(img2>=start)[0],end)]=scale_arg
	return(img2)
	
def gaa(img_gray):
	rawgray=np.float32(img_gray)
	imhist,bins = np.histogram(img_gray.flatten(),255,normed=True)
	
	if np.where(imhist[251:255]==0)[0].size >0:
		maxpos=np.where(imhist[251:255]==0)[0][0]+251
	else:
		maxpos=255
		
	if np.where(imhist[0:50]==0)[0].size >0:
		minpos=np.where(imhist[0:50]==0)[0][-1]
	else:
		minpos=0
	rawgray[np.where(rawgray<=bins[minpos])]=0
	rawgray[np.where(rawgray>=bins[np.uint8((3*maxpos+255)/4)])]=bins[-1]
	rawgray1=np.uint8(histeq2(rawgray-rawgray.min(),255))
	rawgray=np.uint8(histeq(rawgray1,255))
	gam2= exposure.adjust_gamma(rawgray1, 1)
	return(gam2)

def image_sep_forlyar_0(img_gray,channel):
	tmparray=np.zeros(img_gray.size*3)
	tmpfig=tmparray.reshape(img_gray.shape[0],img_gray.shape[1],3)
	#mi=raw.min()
	#ma=rawgray.max()
	rawgray=np.float32(img_gray)
	imhist,bins = np.histogram(img_gray.flatten(),255,normed=True)
	
	if np.where(imhist[251:255]==0)[0].size >0:
		maxpos=np.where(imhist[251:255]==0)[0][0]+251
	else:
		maxpos=255
		
	if np.where(imhist[0:50]==0)[0].size >0:
		minpos=np.where(imhist[0:50]==0)[0][-1]
	else:
		minpos=0

	rawgray[np.where(rawgray<=bins[minpos])]=0
	rawgray[np.where(rawgray>=bins[np.uint8((3*maxpos+255)/4)])]=bins[-1]
	rawgray1=np.uint8(histeq2(rawgray-rawgray.min(),255))
	rawgray=np.uint8(histeq(rawgray1,255))
	gam2= exposure.adjust_gamma(rawgray1, 1.3)
	tmpfig[:,:,channel]=np.uint8(gam2)
	#fig_c,purple_c,raw_c,pos,count=image_circle(img_gray, tmpfig, abs(2-channel),253.5)
	raw_c=tmpfig.copy()
	pos=[]
	count=0
	r_radius=[]
	gray2=np.uint8(histeq2(img_gray,255))
	sobelx8u = cv2.Sobel(gray2,cv2.CV_8U,1,0,ksize=3)
	gradient=sobelx8u
	return tmpfig,raw_c,pos,count,rawgray,gradient,r_radius
	
### subfunctions for 2 marked proteins images
def image_sep(img_gray,channel,minsize,maxsize):
	tmparray=np.zeros(img_gray.size*3)
	tmpfig=tmparray.reshape(img_gray.shape[0],img_gray.shape[1],3)
	#mi=raw.min()
	#ma=rawgray.max()
	rawgray=np.float32(img_gray)
	clahe=cv2.createCLAHE(clipLimit=2.0,tileGridSize=(8,8))
	#dst=clahe.apply(np.uint8(img_gray))
	imhist,bins = np.histogram(img_gray.flatten(),255,normed=True)
	
	if np.where(imhist[251:255]==0)[0].size >0:
		maxpos=np.where(imhist[251:255]==0)[0][0]+251
	else:
		maxpos=255
		
	if np.where(imhist[0:50]==0)[0].size >0:
		minpos=np.where(imhist[0:50]==0)[0][-1]
	else:
		minpos=0

	rawgray[np.where(rawgray<=bins[minpos])]=0
	rawgray[np.where(rawgray>=bins[np.uint8((3*maxpos+255)/4)])]=bins[-1]
	rawgray1=np.uint8(histeq2(rawgray-rawgray.min(),255))
	rawgray=np.uint8(histeq(rawgray1,255))
	gam2= exposure.adjust_gamma(rawgray1, 1.5)
	dst=clahe.apply(gam2)
	tmpfig[:,:,channel]=np.uint8(dst)
	#fig_c,purple_c,raw_c,pos,count=image_circle(img_gray, tmpfig, abs(2-channel),253.5)
	raw_c,pos,count,gradient,r_radius,gray_c=image_circle(dst, tmpfig, channel,253.5,minsize,maxsize)
	return tmpfig,raw_c,pos,count,dst,gradient,r_radius,gray_c
 
### subfunctions for 2 marked proteins images
def image_sep3(img_gray,channel):
	tmparray=np.zeros(img_gray.size*3)
	tmpfig=tmparray.reshape(img_gray.shape[0],img_gray.shape[1],3)
	#mi=raw.min()
	#ma=rawgray.max()
	rawgray=np.float32(img_gray)
	imhist,bins = np.histogram(img_gray.flatten(),255,normed=True)
	
	if np.where(imhist[241:255]==0)[0].size >0:
		maxpos=np.where(imhist[241:255]==0)[0][0]+241
	else:
		maxpos=255
		
	if np.where(imhist[0:50]==0)[0].size >0:
		minpos=np.where(imhist[0:50]==0)[0][-1]
	else:
		minpos=0

	rawgray[np.where(rawgray<=bins[minpos])]=0
	rawgray[np.where(rawgray>=bins[np.uint8((3*maxpos+255)/4)])]=bins[-1]
	rawgray1=np.uint8(histeq2(rawgray-rawgray.min(),255))
	rawgray=np.uint8(histeq2(rawgray1,255))
	gam2= exposure.adjust_gamma(rawgray1, 1.3)
	tmpfig[:,:,channel]=np.uint8(gam2)
	
	return tmpfig
 	

def check_point(pos,set):
	for p in set:
		x1=pos[0]
		x2=p[0]
		y1=pos[1]
		y2=p[1]
		dis=np.linalg.norm(np.array([x1,y1])-np.array([x2,y2]))
		if dis<6:
			return True
			
	return False


def det(gray,point,val):
	h=np.uint(point[0])
	w=np.uint(point[1])
	m=np.uint([h-1,h,h+1])
	n=np.uint([w-1,w,w+1])
	k=0
	for i in m:
		for j in n:
			if(gray[i,j]>val):
				k=k+1
	if k>4:
		return True
	else:
		return False

def rect_angle_sel(grayimg,pos,rr):
	n=grayimg[np.uint(pos[0]-rr):np.uint(pos[0]+rr),np.uint(pos[1]-rr):np.uint(pos[1]+rr)]
	return n
		
def max_near(gray,point):
	h=np.uint(point[0])
	w=np.uint(point[1])
	m=np.uint([h-1,h,h+1])
	n=np.uint([w-1,w,w+1])
	k=0
	pos2=[]
	for i in m:
		for j in n:
			if(gray[i,j]>k):
				k=gray[i,j]
				pos2=[i,j]
	return k,pos2

def max_near2(gray,point):
	h=np.uint(point[0])
	w=np.uint(point[1])
	m=np.uint([h-1,h,h+1])
	n=np.uint([w-1,w,w+1])
	k=0
	pos2=[]
	for i in m:
		for j in n:
			if(gray[i,j]>k):
				k=gray[i,j]
				pos2=[i,j]
	return k	
	

def image_circle(img_gray,tmpfig,channel,cutoff,sizemin,sizemax):
	''' channel 0 is red channel, channel=1 is green channel; cutoff value maybed changed '''
	#img_gray=cv2.imread("img_000000000_YFP_007.tif",-1)
	gray2=gaa(img_gray)
	m=np.zeros(img_gray.size*3)
	m.shape=(img_gray.shape[0],img_gray.shape[1],3)
	gray2=np.uint8(histeq2(img_gray,255))
	m[:,:,channel]=gray2
	#m=np.uint8(m)

	img=m.copy()

	k5=cv2.GaussianBlur(gray2,(7,7),0)
	#(_, thresh) = cv2.threshold(k5, 45,255, cv2.THRESH_BINARY)
	 
	#sobelx8u = cv2.Sobel(k5,cv2.CV_8U,1,0,ksize=5)
	#sobelx8u2 = cv2.Sobel((255-k5),cv2.CV_8U,1,0,ksize=5)

	#save_img(m,"1.tif")

	sobelx8u= cv2.Sobel(k5,cv2.CV_8U,1,0,ksize=3)
	params = cv2.SimpleBlobDetector_Params()
	# Change thresholds
	params.minThreshold = 1
	 
	params.maxThreshold = 254

	# Filter by Area.
	params.filterByArea =True
	params.minArea =sizemin
	params.maxArea=sizemax
	#Filter by Circularity
	params.filterByCircularity = True
	params.minCircularity = 0.001
	# Filter by color

	params.filterByColor = True
	params.blobColor =255
	# Filter by Convexity
	params.filterByConvexity = 0
	params.minConvexity = 0.01
	 
	# Filter by Inertia
	params.filterByInertia = True
	params.minInertiaRatio = 0.05
	# Create a detector with the parameters
	detector = cv2.SimpleBlobDetector_create(params)
	(_, thresh) = cv2.threshold(sobelx8u, 5,1, cv2.THRESH_BINARY)
	val=np.multiply(thresh,np.uint8(sobelx8u))


	keypoints = detector.detect(val)
	len(keypoints)
	k=len(keypoints)
		
	img3=img.copy()
	real_radius=[]
	realpos=[]
	real_keypoint=[]
	for i in range(k):
		c=(np.uint(keypoints[i].pt[1]),np.uint(keypoints[i].pt[0]))
		if max_near(thresh,c)[0]>0 and max_near(sobelx8u,c)[0]>30 and max_near(gray2,c)[0]>30:
			img3=cv2.circle(img3,(c[1],c[0]),15,(0,abs(channel-1)*65535,abs(2-channel)*65535),2)
			real_radius.append(keypoints[i].size)
			realpos.append(c)
			real_keypoint.append(keypoints[i])
		#cv2.imwrite('1.tif',img3)

	len(real_keypoint)

	im_with_keypoints = cv2.drawKeypoints(gray2, real_keypoint, np.array([]), (0,0,255), cv2.DRAW_MATCHES_FLAGS_DRAW_RICH_KEYPOINTS)


	return img3,realpos,k,sobelx8u,real_radius,im_with_keypoints
 

def de_PF_point(np_arr,val):
	c=0
	for i in np.uint(np_arr):
		if i>=val:
			c=c+1
	if c>1:
		return True
	else:
		return False
def de_PF_points(np_arr,val):
	c=0
	for i in np.uint(np_arr):
		if i>=val:
			c=c+12
	if c>13:
		return False
	else:
		return True
 
def mark_second_circle(pos,graydata,channel):
	m=np.zeros(graydata.size*3)
	m.shape=(graydata.shape[0],graydata.shape[1],3)
	gray=np.uint8(histeq2(graydata,255))
	m[:,:,channel]=gray
	img=m.copy()
	for i in range(len(pos)):
		c=pos[i]
		img=cv2.circle(img,(np.uint(c[1]),np.uint(c[0])),15,(0,abs(channel-1)*65535,abs(2-channel)*65535),2) 
	return img
	
def mufun(gray,pos):
	k=0
	j=0
	posy=pos[0]
	posx=pos[1]
	while gray[np.uint(posy+k),np.uint(posx+j)]>200:
		k=k+1
		j=j+1
	return(np.uint(posy+k),np.uint(posx+j))

def match_block(tar_gray,key_gray):
	w=key_gray.shape[1]
	h=key_gray.shape[0]
	for i in range(1+tar_gray.shape[0]-h):
		for j in range(1+tar_gray.shape[1]-w):
			bloc=tar_gray[i:(i+h),j:(j+w)]
			if sum(sum(bloc==key_gray))==np.uint(w*h):
				return True,[i+h/2,j+w/2]
	return False,[i+h/2,j+w/2]


def det2(gray,point):
	h=np.uint(point[0])
	w=np.uint(point[1])
	tar_blo1=np.uint(np.array([[255,255,255],[255,255,255],[255,255,255]]))
	tar_blo2=np.uint(np.array([[0,0,0],[0,0,0],[0,0,0]]))
	#tar_blo1=np.uint(np.array([[255,255],[255,255]]))
	#tar_blo2=np.uint(np.array([[0,0],[0,0]]))
	for i in range(5):
		for l in range(3):
			for m in range(3):
				tar_blo1[l,m]=255-i
				a,b=match_block(gray[np.uint(h-5):np.uint(h),np.uint(w-5):np.uint(w)],tar_blo1)
				c,d=match_block(gray[np.uint(h):np.uint(h+5),np.uint(w):np.uint(w+5)],tar_blo2)
				if a and c:
					return True
	else:
		return False
		
def det_fir(gray,point):
	h=np.uint(point[0])
	w=np.uint(point[1])
	#tar_blo1=np.uint(np.array([[255,255,255],[255,255,255],[255,255,255]]))
	#tar_blo2=np.uint(np.array([[0,0,0],[0,0,0],[0,0,0]]))
	tar_blo1=np.uint(np.array([[255,255],[255,255]]))
	tar_blo2=np.uint(np.array([[0,0],[0,0]]))
	if match_block(gray[np.uint(h-5):np.uint(h),np.uint(w-5):np.uint(w)],tar_blo1) and match_block(gray[np.uint(h):np.uint(h+5),np.uint(w):np.uint(w+5)],tar_blo2):
		return True
	else:
		return False
