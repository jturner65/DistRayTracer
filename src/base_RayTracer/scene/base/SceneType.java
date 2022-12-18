package base_RayTracer.scene.base;

import java.util.HashMap;
import java.util.Map;

/**
 * This enum describes the type of scene to be created
 * @author John Turner
 *
 */
public enum SceneType {
	Fisheye(0),
	FOVScene(1),
	DepthofFieldScene(2),
	Orthographic(3);
	private int value;
	
	private String[] _typeExplanation = new String[] {
			"Fisheye Lense Scene",
			"Field of View Scene",
			"FOV Scene with Depth Of Field",
			"Orthographic Scene"};
	private static String[] _typeName = new String[] {
		"Fisheye","FOV","FOV_Depth","Ortho"	
	};
	
	public static String[] getListOfTypes() {return _typeName;}
	private static Map<Integer, SceneType> map = new HashMap<Integer, SceneType>(); 
	static { for (SceneType enumV : SceneType.values()) { map.put(enumV.value, enumV);}}
	private SceneType(int _val){value = _val;} 
	public int getVal(){return value;}
	public static SceneType getVal(int idx){return map.get(idx);}
	public static int getNumVals(){return map.size();}						//get # of values in enum
	public String getName() {return _typeName[value];}
	@Override
    public String toString() { return ""+_typeExplanation[value] + "("+value+")"; }	
    public String toStrBrf() { return ""+_typeExplanation[value]; }			
	
}//enum SceneType
