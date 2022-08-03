package uamds.vis;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;
import java.util.function.BiConsumer;
import java.util.function.DoubleBinaryOperator;
import java.util.stream.Collectors;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;

import hageldave.jplotter.canvas.BlankCanvas;
import hageldave.jplotter.canvas.BlankCanvasFallback;
import hageldave.jplotter.canvas.JPlotterCanvas;
import hageldave.jplotter.interaction.kml.CoordSysPanning;
import hageldave.jplotter.interaction.kml.CoordSysScrollZoom;
import hageldave.jplotter.interaction.kml.KeyMaskListener;
import hageldave.jplotter.misc.Contours;
import hageldave.jplotter.renderables.Lines;
import hageldave.jplotter.renderables.Lines.SegmentDetails;
import hageldave.jplotter.renderables.Triangles;
import hageldave.jplotter.renderables.Triangles.TriangleDetails;
import hageldave.jplotter.renderers.ChainedRenderer;
import hageldave.jplotter.renderers.CompleteRenderer;
import hageldave.jplotter.renderers.CoordSysRenderer;
import hageldave.jplotter.renderers.Renderer;
import hageldave.jplotter.util.PickingRegistry;
import uamds.misc.Utils2;
import uamds.other.Utils;

public class IsoLinesPlot {

	public final CompleteRenderer completeRendererBG = new CompleteRenderer();
	public final CompleteRenderer completeRendererFG = new CompleteRenderer();
	public final ChainedRenderer contentRenderer = completeRendererBG.withAppended(completeRendererFG);
	public final CoordSysRenderer coordsys;
	public ArrayList<Lines> isoLines = new ArrayList<>();
	public ArrayList<Triangles> isoBands = new ArrayList<>();
	
	protected PickingRegistry<Integer> pickingRegistry = new PickingRegistry<>();
	protected TreeMap<Integer, Integer> idx2pickcolor = new TreeMap<>();
	protected List<BiConsumer<Integer, MouseEvent>> pickingClickListeners = new LinkedList<>();
	protected List<BiConsumer<Integer, MouseEvent>> pickingMoveListeners = new LinkedList<>();
	
	public final JPlotterCanvas canvas;
	public JFrame frame;
	
	protected int isoLinesSpaceResolution = 64;
	
	public IsoLinesPlot() {
		this(new BlankCanvasFallback(), new CoordSysRenderer());
		canvas.setRenderer(coordsys);
		createMousePickColorClickCapabilities();
	}
	
	public IsoLinesPlot(JPlotterCanvas canvas, CoordSysRenderer csr) {
		this.canvas = canvas;
		this.coordsys = csr;
		Renderer content = coordsys.getContent();
		Renderer ownContent = contentRenderer;
		coordsys.setContent(content == null ? ownContent : content.withAppended(ownContent));
	}
	
	public int getNumIsoLineInstances() {
		return isoLines.size();
	}
	
	public JPlotterCanvas getCanvas() {
		return canvas;
	}
	
	public int addEmpty() {
		Lines lines = new Lines();
		isoLines.add(lines);
		int idx = isoLines.size()-1;
		
		completeRendererFG.addItemToRender(lines);
		int pickcolor = pickingRegistry.register(idx);
		idx2pickcolor.put(idx, pickcolor);
		
		Triangles tris = new Triangles();
		isoBands.add(tris);
		completeRendererFG.addItemToRender(tris);
		
		return idx;
	}
	
	public int addIsoLines(DoubleBinaryOperator fx, Rectangle2D view, double[] levels, int[] colors) {
		int idx = addEmpty();
		updateIsoLines(idx, fx, view, levels, colors);
		return idx;
	}
	
	public void purgeLines(int idx) {
		SwingUtilities.invokeLater(()->{
			isoLines.get(idx).removeAllSegments();
		});
		if(canvas != null)
			canvas.scheduleRepaint();
	}
	
	public void updateIsoLines(int idx, DoubleBinaryOperator fx,  Rectangle2D view, double[] levels, int[] colors) {
		Lines lines = isoLines.get(idx);
		Triangles bands = isoBands.get(idx);
		double xmin= view.getMinX();
		double xmax= view.getMaxX();
		double ymin= view.getMinY();
		double ymax= view.getMaxY();
		int resolution = isoLinesSpaceResolution;
		double[][][] XY = Utils2.sampling2DCoordinates(xmin, xmax, ymin, ymax, resolution);
		double[][] Z = Utils2.sample2D(fx, XY[0], XY[1]);
		List<SegmentDetails> segments = new LinkedList<>();
		List<TriangleDetails> bandtris = new LinkedList<>();
		for(int i=0; i<levels.length; i++) {
			List<SegmentDetails> segs = Contours.computeContourLines(XY[0], XY[1], Z, levels[i], colors[i]);
			segments.addAll(segs);
			
//			List<TriangleDetails> tris = Contours.computeContourBands(XY[0], XY[1], Z, i==0 ? 1000:levels[i-1], levels[i], colors[i], colors[i]);
//			bandtris.addAll(tris);
		}
		hageldave.jplotter.util.Utils.execOnAWTEventDispatch(()->{
			lines.removeAllSegments().getSegments().addAll(segments);
			int pickColor = idx2pickcolor.get(idx);
			lines.getSegments().forEach(seg->seg.setPickColor(pickColor));
			
			bands.removeAllTriangles().getTriangleDetails().addAll(bandtris);
		});
		if(canvas != null)
			canvas.scheduleRepaint();
	}
	
	protected void createMousePickColorClickCapabilities() {
		MouseAdapter ma = new MouseAdapter() {
			@Override
			public void mouseClicked(MouseEvent e) {
				int pixel = canvas.getPixel(e.getX(), e.getY(), true, 5);
				Integer idx = pickingRegistry.lookup(pixel);
				notifyPickingClick(idx, e);
			}
			
			@Override
			public void mouseMoved(MouseEvent e) {
				int pixel = canvas.getPixel(e.getX(), e.getY(), true, 5);
				Integer idx = pickingRegistry.lookup(pixel);
				notifyPickingMove(idx, e);
			}
		};
		this.canvas.asComponent().addMouseListener(ma);
		this.canvas.asComponent().addMouseMotionListener(ma);
	}
	
	protected void notifyPickingClick(Integer idx, MouseEvent e) {
		for(var l : pickingClickListeners) {
			l.accept(idx, e);
		}
	}
	
	protected void notifyPickingMove(Integer idx, MouseEvent e) {
		for(var l : pickingMoveListeners) {
			l.accept(idx, e);
		}
	}
	
	public void addPickingClickListener(BiConsumer<Integer, MouseEvent> listener) {
		this.pickingClickListeners.add(listener);
	}
	
	public void addPickingMoveListener(BiConsumer<Integer, MouseEvent> listener) {
		this.pickingMoveListeners.add(listener);
	}
	
	public JFrame display(String title) {
		frame = new JFrame(title);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.setContentPane(new JPanel(new BorderLayout()));
		
		canvas.addCleanupOnWindowClosingListener(frame);
		
		frame.getContentPane().add(canvas.asComponent(), BorderLayout.CENTER);
		canvas.asComponent().setPreferredSize(new Dimension(480, 480));
		canvas.asComponent().setBackground(Color.white);

		new CoordSysScrollZoom(canvas, coordsys, new KeyMaskListener(0))
		.setZoomFactor(1.25)
		.register();
		new CoordSysPanning(canvas, coordsys, new KeyMaskListener(0)).register();
		
		SwingUtilities.invokeLater(()->{
			frame.pack();
			frame.setVisible(true);
		});
		
		return frame;
		
	}
	
	public synchronized void highLight(int... idx) {
		for(Lines l : isoLines) {
			completeRendererBG.lines.removeItemToRender(l);
			completeRendererFG.lines.removeItemToRender(l);
		}
		Set<Integer> toHighLight = Arrays.stream(idx).boxed().collect(Collectors.toSet());
		if(toHighLight.isEmpty()) {
			// reset
			for(Lines l : isoLines) {
				completeRendererFG.addItemToRender(l);
				l.setGlobalSaturationMultiplier(1.0);
				l.setGlobalAlphaMultiplier(1.0);
			}
		} else {
			for(int i = 0; i<isoLines.size(); i++) {
				Lines l = isoLines.get(i);
				if(toHighLight.contains(i)) {
					completeRendererFG.addItemToRender(l);
					l.setGlobalSaturationMultiplier(1.0);
					l.setGlobalAlphaMultiplier(1.0);
				} else {
					completeRendererBG.addItemToRender(l);
					l.setGlobalSaturationMultiplier(0.2);
					l.setGlobalAlphaMultiplier(0.2);
				}
			}
		}
		canvas.scheduleRepaint();
	}

	public synchronized void emphasize(int... idx) {
		for(Lines l : isoLines) {
			l.setGlobalThicknessMultiplier(1);
		}
		for(int i : idx) {
			isoLines.get(i).setGlobalThicknessMultiplier(2);
		}
		canvas.scheduleRepaint();
	}
	
}
