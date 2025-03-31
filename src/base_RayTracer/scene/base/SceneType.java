package base_RayTracer.scene.base;

import java.util.HashMap;
import java.util.Map;

/**
 * This enum describes the type of scene to be created
 * @author John Turner
 *
 */
public enum SceneType {
	Fisheye,
	FOVScene,
	DepthofFieldScene,
	Orthographic;
	
	private final String[] _typeExplanation = new String[] {
			"Fisheye Lense Scene",
			"Field of View Scene",
			"FOV Scene with Depth Of Field",
			"Orthographic Scene"};
	private static final String[] _typeName = new String[] {
		"Fisheye","FOV","FOV_Depth","Ortho"	
	};
	
	public static String[] getListOfTypes() {return _typeName;}
	private static Map<Integer, SceneType> map = new HashMap<Integer, SceneType>(); 
	static { for (SceneType enumV : SceneType.values()) { map.put(enumV.ordinal(), enumV);}}
	public int getVal(){return ordinal();}
	public static SceneType getEnumByIndex(int idx){return map.get(idx);}
	public static SceneType getEnumFromValue(int idx){return map.get(idx);}
	public static int getNumVals(){return map.size();}						//get # of values in enum
	public String getName() {return _typeName[ordinal()];}
	@Override
    public String toString() { return ""+_typeExplanation[ordinal()] + "("+ordinal()+")"; }	
    public String toStrBrf() { return ""+_typeExplanation[ordinal()]; }			
	
}//enum SceneType
