package base_RayTracer.ui;

import java.io.File;
import java.util.ArrayList;
import java.util.TreeMap;

import base_Render_Interface.IRenderInterface;
import base_RayTracer.scene.base.Base_Scene;
import base_RayTracer.ui.base.Base_RayTracerUIUpdater;
import base_RayTracer.ui.base.Base_RayTracerWin;
import base_UI_Objects.GUI_AppManager;

/**
 * class to hold 2-D ray tracer experiment - bunch of circles, shoot rays and plot their traversal
 * @author john
 *
 */
public class RayTracer2DWin extends Base_RayTracerWin {
	
	/////////////
	// ui objects 
	////////////

		
	public RayTracer2DWin(IRenderInterface _p, GUI_AppManager _AppMgr, int _winIdx) {
		super(_p, _AppMgr, _winIdx);
		super.initThisWin(false);
		setNewUIValue(gIDX_SceneCols, 600);setUIWinVals(gIDX_SceneCols);
		setNewUIValue(gIDX_SceneRows, 600);setUIWinVals(gIDX_SceneRows);
		startRayTrace();
	}//ctor
			
	@Override
	protected void initMe_Indiv() {}//
	
	//initialize all UI buttons here
	@Override
	public int initAllPrivBtns_Indiv(ArrayList<Object[]> tmpBtnNamesArray) {
		return numPrivFlags;	
	}

	@Override
	protected void setupGUIObjsAras_Indiv(TreeMap<Integer, Object[]> tmpUIObjArray,
			TreeMap<Integer, String[]> tmpListObjVals) {
	}

	@Override
	protected void setPrivFlags_Indiv(int idx, boolean val) {		
	}
	/**
	 * This function would provide an instance of the override class for base_UpdateFromUIData, which would
	 * be used to communicate changes in UI settings directly to the value consumers.
	 */
	@Override
	protected Base_RayTracerUIUpdater buildUIDataUpdateObject_Indiv() {
		return new Base_RayTracerUIUpdater(this);
	}
	/**
	 * This function is called on ui value update, to pass new ui values on to window-owned consumers
	 */
	@Override
	protected final void updateCalcObjUIVals() {}

	/**
	 * UI code-level Debug mode functionality. Called only from flags structure
	 * @param val
	 */
	@Override
	protected final void handleDispFlagsDebugMode_Indiv(boolean val) {}
	
	/**
	 * Application-specific Debug mode functionality (application-specific). Called only from privflags structure
	 * @param val
	 */
	@Override
	protected final void handlePrivFlagsDebugMode_Indiv(boolean val) {	}

	@Override
	protected boolean setUI_IntValsCustom_Indiv(int UIidx, int ival, int oldVal) {	return false;}

	@Override
	protected boolean setUI_FloatValsCustom_Indiv(int UIidx, float val, float oldVal) {	return false;}

	
	@Override
	protected int[] getFlagIDXsToInitToTrue_Indiv(int[] flagsToInit) {
		return flagsToInit;
	}

	//check whether the mouse is over a legitimate map location
	@Override
	protected boolean chkMouseClick2D(int mouseX, int mouseY, int btn){		
		return false;
	}//chkMouseOvr
	
	//check whether the mouse is over a legitimate map location
	@Override
	protected boolean chkMouseMoveDragState2D(int mouseX, int mouseY, int btn){		
		return false;
	}//chkMouseOvr
	
	//check whether the mouse is over a legitimate map location
	@Override
	protected void setMouseReleaseState2D(){}//chkMouseOvr

	@Override
	protected void drawMe_Indiv(float animTimeMod) {
		ri.pushMatState();
		float[] loc = getLocUpperCrnr();
		ri.translate(loc[0],loc[1],0);
		ri.pushMatState();
		ri.disableLights();
		Base_Scene s = getCurrScene();
		if(s!=null) {s.draw();}	
		ri.popMatState();
		ri.popMatState();
	}

	@Override
	protected void drawRightSideInfoBarPriv_Indiv(float modAmtMillis) {		
	}

	@Override
	protected void drawCustMenuObjs_Indiv() {

	}

	@Override
	protected void setVisScreenDimsPriv_Indiv() {}

	@Override
	public void hndlFileLoad_Indiv(File file, String[] vals, int[] stIdx) {}
	
	
}//RayTracer2DWin
