package uamds.vis;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;

import hageldave.jplotter.canvas.BlankCanvas;
import hageldave.jplotter.canvas.BlankCanvasFallback;
import hageldave.jplotter.canvas.JPlotterCanvas;
import hageldave.jplotter.interaction.kml.CoordSysViewSelector;
import hageldave.jplotter.interaction.kml.CoordSysPanning;
import hageldave.jplotter.interaction.kml.CoordSysScrollZoom;
import hageldave.jplotter.interaction.kml.KeyMaskListener;
import hageldave.jplotter.renderers.CoordSysRenderer;

public class CoordSysDisplay {

	public CoordSysRenderer coordsys = new CoordSysRenderer();
	public JPlotterCanvas canvas;
	public JFrame frame = new JFrame();
	
	
	public CoordSysDisplay(boolean opengl) {
		canvas = opengl ? new BlankCanvas():new BlankCanvasFallback();
		canvas.addCleanupOnWindowClosingListener(frame);
		canvas.setRenderer(coordsys);
		frame.setContentPane(new JPanel(new BorderLayout()));
		frame.getContentPane().add(canvas.asComponent(), BorderLayout.CENTER);
		canvas.asComponent().setPreferredSize(new Dimension(480, 480));
		canvas.asComponent().setBackground(Color.white);
		new CoordSysScrollZoom(canvas, coordsys, new KeyMaskListener(0)).register();
		new CoordSysPanning(canvas, coordsys).register();
		new CoordSysViewSelector(canvas, coordsys) {
			@Override
			public void areaSelected(double minx, double miny, double maxx, double maxy) {
				coordsys.setCoordinateView(minx, miny, maxx, maxy);
			}
		}.register();
	}
	
	public JFrame setupAndShow() {
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		SwingUtilities.invokeLater(()->{
			frame.pack();
			frame.setVisible(true);
		});
		
		return frame;
	}
}
