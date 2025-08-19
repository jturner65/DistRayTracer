package base_RayTracer.ui;

import java.io.File;
import java.util.LinkedHashMap;

import base_RayTracer.scene.base.Base_Scene;
import base_RayTracer.ui.base.Base_RayTracerUIUpdater;
import base_RayTracer.ui.base.Base_RayTracerWin;
import base_Render_Interface.IGraphicsAppInterface;
import base_UI_Objects.GUI_AppManager;
import base_UI_Objects.windowUI.base.GUI_AppWinVals;
import base_UI_Objects.windowUI.uiObjs.base.GUIObj_Params;

/**
 * class to hold 2-D ray tracer experiment - bunch of circles, shoot rays and plot their traversal
 * @author john
 *
 */
public class RayTracer2DWin extends Base_RayTracerWin {
    
    /////////////
    // ui objects 
    ////////////

        
    public RayTracer2DWin(IGraphicsAppInterface _p, GUI_AppManager _AppMgr, int _winIdx) {
        super(_p, _AppMgr, _winIdx);
        
    }//ctor
            
    @Override
    protected void initMe_Indiv() {        
        uiMgr.forceNewUIValue(gIDX_SceneCols, 600);uiMgr.updateOwnerWithUIVal(gIDX_SceneCols);
        uiMgr.forceNewUIValue(gIDX_SceneRows, 600);uiMgr.updateOwnerWithUIVal(gIDX_SceneRows);
        startRayTrace(AppMgr.isGlblDebug());
    }//
    /**
     * This function implements the instantiation of a child window owned by this window, if such exists.
     * The implementation should be similar to how the main windows are implemented in GUI_AppManager::initAllDispWindows.
     * If no child window exists, this implementation of this function can be empty
     * If a child window is instantiated, it MUST have its init called (childWin.initThisWin(false))
     * 
     * @param GUI_AppWinVals the window control values for the child window.
     */
    @Override
    protected final void buildAndSetChildWindow_Indiv(GUI_AppWinVals _appVals) {}  
    
    /**
     * Retrieve the total number of defined privFlags booleans (application-specific state bools and interactive buttons)
     */
    @Override
    public int getTotalNumOfPrivBools() { return numPrivFlags;}

    /**
     * Build all UI objects to be shown in left side bar menu for this window. This is the first child class function called by initThisWin
     * @param tmpUIObjMap : map of GUIObj_Params, keyed by unique string, with values describing the UI object
     *             - The object IDX                   
     *          - A double array of min/max/mod values                                                   
     *          - The starting value                                                                      
     *          - The label for object                                                                       
     *          - The object type (GUIObj_Type enum)
     *          - A boolean array of behavior configuration values : (unspecified values default to false)
     *               idx 0: value is sent to owning window,  
     *               idx 1: value is sent on any modifications (while being modified, not just on release), 
     *               idx 2: changes to value must be explicitly sent to consumer (are not automatically sent),
     *          - A boolean array of renderer format values :(unspecified values default to false) - Behavior Boolean array must also be provided!
     *                 - Should be multiline
     *                 - One object per row in UI space (i.e. default for multi-line and btn objects is false, single line non-buttons is true)
     *                 - Force this object to be on a new row/line (For side-by-side layouts)
     *                 - Text should be centered (default is false)
     *                 - Object should be rendered with outline (default for btns is true, for non-buttons is false)
     *                 - Should have ornament
     *                 - Ornament color should match label color 
     */
    @Override
    protected final void setupGUIObjsAras_Indiv(LinkedHashMap<String, GUIObj_Params> tmpUIObjMap) {}

    /**
     * Build all UI buttons to be shown in left side bar menu for this window. This is for instancing windows to add to button region
     * @param firstIdx : the first index to use in the map/as the objIdx
     * @param tmpUIBoolSwitchObjMap : map of GUIObj_Params to be built containing all flag-backed boolean switch definitions, keyed by sequential value == objId
     *                 the first element is the object index
     *                 the second element is true label
     *                 the third element is false label
     *                 the final element is integer flag idx 
     */
    @Override
    protected final void setupGUIBoolSwitchAras_Indiv(int firstIdx, LinkedHashMap<String, GUIObj_Params> tmpUIBoolSwitchObjMap){}

    @Override
    protected void setPrivFlags_Indiv(int idx, boolean val) {}
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
    protected final void handlePrivFlagsDebugMode_Indiv(boolean val) {    }

    @Override
    protected boolean setUI_IntValsCustom_Indiv(int UIidx, int ival, int oldVal) {    return false;}

    @Override
    protected boolean setUI_FloatValsCustom_Indiv(int UIidx, float val, float oldVal) {    return false;}

    
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
    protected void drawMe_Indiv(float animTimeMod, boolean isGlblAppDebug) {
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
    protected void drawRightSideInfoBarPriv_Indiv(float modAmtMillis, boolean isGlblAppDebug) {        
    }

    @Override
    protected void drawCustMenuObjs_Indiv(boolean isGlblAppDebug) {

    }

    @Override
    protected void setVisScreenDimsPriv_Indiv() {}

    @Override
    public void hndlFileLoad_Indiv(File file, String[] vals, int[] stIdx) {}
    
    
}//RayTracer2DWin
