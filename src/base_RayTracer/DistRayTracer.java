package base_RayTracer;

import java.util.HashMap;

import base_RayTracer.ui.RayTracer2DWin;
import base_UI_Objects.GUI_AppManager;
import base_UI_Objects.windowUI.sidebar.SidebarMenu;
import base_Utils_Objects.io.messaging.MsgCodes;

/**
 * this file can act as a stub for the ray tracer and can launch it
 * @author John Turner
 *
 */
public class DistRayTracer extends GUI_AppManager {
	//project-specific variables
	public final String prjNmLong = "Testbed for base ray tracer";
	public final String prjNmShrt = "RayTracerBaseExp";
	public final String projDesc = "Testbed for distribution base ray tracer";


	/**
	 * idx's in dispWinFrames for each window - 0 is always left side menu window
	 * Side menu is dispMenuIDX == 0
	 */
	private static final int disp2DRayTracerIDX = 1;
	/**
	 * # of visible windows including side menu (always at least 1 for side menu)
	 */
	private static final int numVisWins = 2;
	//don't use sphere background for this program
	private boolean useSphereBKGnd = false;	
	
	private String bkSkyBox = "bkgrndTex.jpg";

	private final int[] bground = new int[]{244,244,244,255};		//bground color

		
	public static void main(String[] passedArgs) {		
		DistRayTracer me = new DistRayTracer();
		DistRayTracer.invokeProcessingMain(me, passedArgs);		    
	}//main
	
	protected DistRayTracer(){super();}

	@Override
	protected boolean showMachineData() {return true;}
	/**
	 * Set various relevant runtime arguments in argsMap
	 * @param _passedArgs command-line arguments
	 */
	@Override
	protected HashMap<String,Object> setRuntimeArgsVals(HashMap<String, Object> _passedArgsMap) {
		return  _passedArgsMap;
	}
	
	/**
	 * Called in pre-draw initial setup, before first init
	 * potentially override setup variables on per-project basis.
	 * Do not use for setting background color or Skybox anymore.
	 *  	(Current settings in my_procApplet) 	
	 *  	strokeCap(PROJECT);
	 *  	textSize(txtSz);
	 *  	textureMode(NORMAL);			
	 *  	rectMode(CORNER);	
	 *  	sphereDetail(4);	 * 
	 */
	@Override
	protected void setupAppDims_Indiv() {}
	@Override
	protected boolean getUseSkyboxBKGnd(int winIdx) {	return useSphereBKGnd;}
	@Override
	protected String getSkyboxFilename(int winIdx) {	return bkSkyBox;}
	@Override
	protected int[] getBackgroundColor(int winIdx) {return bground;}
	@Override
	protected int getNumDispWindows() {	return numVisWins;	}
	
	/**
	 * whether or not we want to restrict window size on widescreen monitors
	 * 
	 * @return 0 - use monitor size regardless
	 * 			1 - use smaller dim to be determine window 
	 * 			2+ - TBD
	 */
	@Override
	protected int setAppWindowDimRestrictions() {	return 1;}		
	
	@Override
	public String getPrjNmShrt() {		return prjNmShrt;}
	@Override
	public String getPrjNmLong() {		return prjNmLong;}
	@Override
	public String getPrjDescr() {		return projDesc;}

	/**
	 * Set minimum level of message object console messages to display for this application. If null then all messages displayed
	 * @return
	 */
	@Override
	protected final MsgCodes getMinConsoleMsgCodes() {return MsgCodes.info1;}
	/**
	 * Set minimum level of message object log messages to save to log for this application. If null then all messages saved to log.
	 * @return
	 */
	@Override
	protected final MsgCodes getMinLogMsgCodes() {return null;}
	
	/**
	 * determine which main flags to show at upper left of menu 
	 */
	@Override
	protected void initBaseFlags_Indiv() {
		setBaseFlagToShow_debugMode(true);
		setBaseFlagToShow_saveAnim(true); 
		setBaseFlagToShow_runSim(false);
		setBaseFlagToShow_singleStep(false);
		setBaseFlagToShow_showRtSideMenu(true);
		setBaseFlagToShow_showDrawableCanvas(false);
	}

	
	@Override
	protected void initAllDispWindows() {
		showInfo = true;
		//titles and descs, need to be set before sidebar menu is defined
		String[] _winTitles = new String[]{"","2D Ray Tracer"},
				_winDescr = new String[] {"", "2D ray tracing renderer."};

		//instanced window dims when open and closed - only showing 1 open at a time - and init cam vals
		float[][] _floatDims  = getDefaultWinAndCameraDims();	

		//application-wide menu button bar titles and button names
		String[] menuBtnTitles = new String[]{"Special Functions 1","Special Functions 2"};
		String[][] menuBtnNames = new String[][] { // each must have literals for every button defined in side bar menu, or ignored
			{"Func 1", "Func 2","Func 3"},	//row 1
			{"Func 1", "Func 2", "Func 3", "Func 4"}};	//row 1
		String[] dbgBtnNames = new String[] {"Debug 0","Debug 1","Debug 2","Debug 3","Debug 4"};
		buildSideBarMenu(_winTitles, menuBtnTitles, menuBtnNames, dbgBtnNames, false, true);		
		//define windows
		/**
		 *  _winIdx The index in the various window-descriptor arrays for the dispWindow being set
		 *  _title string title of this window
		 *  _descr string description of this window
		 *  _dispFlags Essential flags describing the nature of the dispWindow for idxs : 
		 * 		0 : dispWinIs3d, 
		 * 		1 : canDrawInWin; 
		 * 		2 : canShow3dbox (only supported for 3D); 
		 * 		3 : canMoveView
		 *  _floatVals an array holding float arrays for 
		 * 				rectDimOpen(idx 0),
		 * 				rectDimClosed(idx 1),
		 * 				initCameraVals(idx 2)
		 *  _intClrVals and array holding int arrays for
		 * 				winFillClr (idx 0),
		 * 				winStrkClr (idx 1),
		 * 				winTrajFillClr(idx 2),
		 * 				winTrajStrkClr(idx 3),
		 * 				rtSideFillClr(idx 4),
		 * 				rtSideStrkClr(idx 5)
		 *  _sceneCenterVal center of scene, for drawing objects (optional)
		 *  _initSceneFocusVal initial focus target for camera (optional)
		 */
		
		//ray tracer window
		int wIdx = disp2DRayTracerIDX;
		setInitDispWinVals(wIdx, _winTitles[wIdx], _winDescr[wIdx], new boolean[]{false,false,false,false}, _floatDims,
				new int[][] {new int[]{20,30,10,255}, new int[]{255,255,255,255},
					new int[]{180,180,180,255}, new int[]{100,100,100,255},
					new int[]{0,0,0,200},new int[]{255,255,255,255}});
		dispWinFrames[wIdx] = new RayTracer2DWin(ri, this, wIdx);

		//specify windows that cannot be shown simultaneously here
		initXORWins(
				new int[]{disp2DRayTracerIDX},
				new int[]{disp2DRayTracerIDX});			
	}
	@Override
	protected void initOnce_Indiv() {
		setWinVisFlag(disp2DRayTracerIDX, true);
	}
	@Override
	protected void initProgram_Indiv() {}	
	
	@Override
	/**
	 * Overriding main because handling 2d + 3d windows
	 */	
	public final void drawMePost_Indiv(float modAmtMillis, boolean is3DDraw){}

	/**
	 * handle key pressed
	 * @param keyVal 0-9, with or without shift ((keyCode>=48) && (keyCode <=57))
	 * @param keyCode actual code of key having been pressed
	 */
	protected void handleKeyPress(char key, int keyCode) {
		switch(key) {
		case ' ' : {toggleSimIsRunning(); break;}							//run sim
		case 'f' : {dispWinFrames[curFocusWin].setInitCamView();break;}					//reset camera
		case 'a' :
		case 'A' : {toggleSaveAnim();break;}						//start/stop saving every frame for making into animation
		case 's' :
		case 'S' : {break;}//saveSS(prjNmShrt);break;}//save picture of current image
		default : {return;}
		}//switch		
	}//handleKeyPress

	//////////////////////////////////////////
	/// graphics and base functionality utilities and variables
	//////////////////////////////////////////
	
	/**
	 * return the number of visible window flags for this application
	 * @return
	 */
	@Override
	public int getNumVisFlags() {return numVisWins;}
	@Override
	//address all flag-setting here, so that if any special cases need to be addressed they can be
	protected void setVisFlag_Indiv(int idx, boolean val ){
		switch (idx){
			case disp2DRayTracerIDX		: {setWinFlagsXOR(disp2DRayTracerIDX, val); break;}		
			default : {break;}
		}			
	}
	
	@Override
	//gives multiplier based on whether shift, alt or cntl (or any combo) is pressed
	public double clickValModMult(){return ((altIsPressed() ? .1 : 1.0) * (shiftIsPressed() ? 10.0 : 1.0));}	
	@Override
	public boolean isClickModUIVal() {
		//TODO change this to manage other key settings for situations where multiple simultaneous key presses are not optimal or conventient
		return altIsPressed() || shiftIsPressed();		
	}
	//get the ui rect values of the "master" ui region (another window) -> this is so ui objects of one window can be made, clicked, and shown displaced from those of the parent windwo
	@Override
	public float[] getUIRectVals_Indiv(int idx, float[] menuClickDim){
			//this.pr("In getUIRectVals for idx : " + idx);
		switch(idx){
			case disp2DRayTracerIDX	: {	return menuClickDim;}			
			default :  return menuClickDim;
			}
	}//getUIRectVals
	@Override
	public void handleShowWin(int btn, int val, boolean callFlags){//display specific windows - multi-select/ always on if sel
		if(!callFlags){//called from setflags - only sets button state in UI to avoid infinite loop
			setMenuBtnState(SidebarMenu.btnShowWinIdx,btn, val);
		} else {//called from clicking on buttons in UI
			//val is btn state before transition 
			boolean bVal = (val == 1?  false : true);
			//each entry in this array should correspond to a clickable window
			setWinVisFlag(winFlagsXOR[btn], bVal);
		}
	}//handleShowWin
	
	/**
	 * any instancing-class-specific colors - colorVal set to be higher than IRenderInterface.gui_OffWhite
	 * @param colorVal
	 * @param alpha
	 * @return
	 */
	@Override
	public int[] getClr_Custom(int colorVal, int alpha) {
		// TODO Auto-generated method stub
		return new int[] {255,255,255,alpha};
	}

	@Override
	public String[] getMouseOverSelBtnLabels() {
		return new String[0];
	}

	@Override
	protected void setSmoothing() {		ri.setSmoothing(0);		}


}//DistRayTracer class
