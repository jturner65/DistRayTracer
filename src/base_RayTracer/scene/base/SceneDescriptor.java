package base_RayTracer.scene.base;

/**
 * This struct-type object will hold all the pertinent values describing a scene.  It is built by the rtxFilleReader
 * 
 * TODO
 * @author John Turner
 *
 */
public class SceneDescriptor {
	/**
	 * Type of scene being built
	 */
	private SceneType sceneType;
	
		
	public SceneDescriptor() {
		// TODO Auto-generated constructor stub
	}
	
	public void setSceneType(SceneType _sceneType) {sceneType = _sceneType;}
	public SceneType getSceneType() {return sceneType;}

}//class SceneDescriptor
