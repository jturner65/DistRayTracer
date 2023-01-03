package base_RayTracer.scene;

import java.util.concurrent.ThreadLocalRandom;

import base_Render_Interface.IRenderInterface;
import base_RayTracer.ray.rayCast;
import base_RayTracer.scene.base.Base_Scene;
import base_RayTracer.ui.base.Base_RayTracerWin;
import base_RayTracer.utils.myRTColor;
import base_Math_Objects.MyMathUtils;
import base_Math_Objects.vectorObjs.doubles.myVector;

public class myOrthoScene extends Base_Scene{
	//width and height of view - for ortho projection. for perspective will be screen width and height
	private double orthoWidth, orthoHeight, orthoScale;			//normalizers for ortho projection
	private double orthoPerRow, orthoPerCol;			//normalizers for ortho projection
	//public List<Future<Boolean>> callOrthoFutures;
	//public List<myOrthoCall> callOrthoCalcs;

	public myOrthoScene(IRenderInterface _p, Base_RayTracerWin _win, String _sceneName,int _numCols, int _numRows, double _orthoWidth, double _orthoHeight) {
		super(_p,_win,_sceneName,_numCols,_numRows);
		setOrthoVals(_orthoWidth,_orthoHeight);	
	}
	
	@Override
	protected final void initVars_Indiv() {}
	
	/**
	 * After image size is changed, recalculate essential scene-specific values that depend on image size
	 */
	@Override
	protected final void setImageSize_Indiv() {
		setOrthoPerPxlVals();		
	}
	
	private void setOrthoVals(double _orthoWidth, double _orthoHeight) {
		orthoWidth = _orthoWidth;
		orthoHeight = _orthoHeight;
		orthoScale = .5 *(orthoWidth + orthoHeight);
		setOrthoPerPxlVals();				
	}
		
	private void setOrthoPerPxlVals() {
		double div = MyMathUtils.min(orthoWidth, orthoHeight) * MyMathUtils.min(sceneCols,sceneRows)/orthoScale;
		orthoPerRow = orthoHeight/div;
		orthoPerCol = orthoWidth/div;		
	}
	
	@Override
	public myRTColor shootMultiRays(double xBseVal, double yBseVal) {
		myRTColor result,aaResultColor;
		double redVal = 0, greenVal = 0, blueVal = 0, rayY, rayX;//,rayYOffset = 1.0/sceneRows, rayXOffset = 1.0/sceneCols;
		//first ray can be straight in
		aaResultColor = reflectRay(new rayCast(this, new myVector(xBseVal,yBseVal,0), new myVector(0,0,-1),0));
		redVal += aaResultColor.x; //(aaResultColor >> 16 & 0xFF)/256.0;//gets red value
		greenVal += aaResultColor.y; // (aaResultColor >> 8 & 0xFF)/256.0;//gets green value
		blueVal += aaResultColor.z;//(aaResultColor & 0xFF)/256.0;//gets blue value	      
		ThreadLocalRandom rand = ThreadLocalRandom.current();
		for(int rayNum = 1; rayNum < numRaysPerPixel; ++rayNum){//vary by +/- .5
			rayY = yBseVal + (orthoPerRow*rand.nextDouble(-.5,.5));
			rayX = xBseVal + (orthoPerCol*rand.nextDouble(-.5,.5));				
			aaResultColor = reflectRay(new rayCast(this, new myVector(rayX,rayY,0), new myVector(0,0,-1),0));
			redVal += aaResultColor.x; //(aaResultColor >> 16 & 0xFF)/256.0;//gets red value
			greenVal += aaResultColor.y; // (aaResultColor >> 8 & 0xFF)/256.0;//gets green value
			blueVal += aaResultColor.z;//(aaResultColor & 0xFF)/256.0;//gets blue value	      
		}//aaliasR
		result = new myRTColor ( redVal/numRaysPerPixel, greenVal/numRaysPerPixel, blueVal/numRaysPerPixel); 
		return result;
	}//shootMultiRays	
	
	@Override
	protected final void renderScene(int stepIter, boolean skipPxl, int[] pixels){
		//index of currently written pixel
		int pixIDX = 0;
		int progressCount = 0;
		myRTColor showColor;
		double rayY, rayX;

		if (numRaysPerPixel == 1){											//single ray into scene per pixel
			for (int row = 0; row < sceneRows; row+=stepIter){
				rayY = orthoPerRow * (-1 * (row - rayYOffset));         
				for (int col = 0; col < sceneCols; col+=stepIter){
					if(skipPxl){skipPxl = false;continue;}			//skip only 0,0 pxl					
					rayX = orthoPerCol * (col - rayXOffset);
					showColor = reflectRay(new rayCast(this,new myVector(rayX,rayY,0), new myVector(0,0,-1),0)); 
					pixIDX = writePxlSpan(showColor.getInt(),row,col,stepIter,pixels);
					progressCount = dispProgressBar(pixIDX, progressCount);//if ((1.0 * pixIDX)/(numPxls) > (progressCount * .02)){System.out.print("-|");++progressCount;}//progressbar   
				}//for col
			}//for row	     
		} else{    //anti aliasing
			for (int row = 0; row < sceneRows; row+=stepIter){
				rayY = orthoPerRow * ((-1 * (row - rayYOffset)) - .5);         
				for (int col = 0; col < sceneCols; col+=stepIter){
					if(skipPxl){skipPxl = false;continue;}			//skip only 0,0 pxl		
					rayX = orthoPerCol * (col - rayXOffset - .5);      
					showColor = shootMultiRays(rayX,rayY); 
					pixIDX = writePxlSpan(showColor.getInt(),row,col,stepIter,pixels);
					progressCount = dispProgressBar(pixIDX, progressCount);//if ((1.0 * pixIDX)/(numPxls) > (progressCount * .02)){System.out.print("-|");++progressCount;}//progressbar  
				}//for col
			}//for row  
		}//if antialiasing			
	}//renderScene	
}//myOrthoScene