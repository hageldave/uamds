package uamds.demo;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;

/**
 * Low fidelity scatter plot component
 */
public class LoFiScatter extends JPanel {
	private static final long serialVersionUID = 1L;
	
	protected ArrayList<double[][]> pointSets = new ArrayList<>();
	protected int pointOpacity = 255;
	protected Rectangle2D contentBounds = calcContentBounds();
	
	public void setPointSets(Collection<double[][]> pointSets) {
		this.pointSets.clear();
		this.pointSets.addAll(pointSets);
		this.contentBounds = calcContentBounds();
		SwingUtilities.invokeLater(this::repaint);
	}
	
	public void setPointOpacity(int pointOpacity) {
		this.pointOpacity = pointOpacity < 0 ? 0:(pointOpacity > 255 ? 255:pointOpacity);
		SwingUtilities.invokeLater(this::repaint);
	}
	
	public Rectangle2D getContentBounds() {
		return contentBounds;
	}

	protected int getColor(int index) {
		int[] colormap = {0xff_66c2a5,
				0xff_fc8d62,
				0xff_8da0cb,
				0xff_e78ac3,
				0xff_a6d854,
				0xff_ffd92f,
				0xff_e5c494,
				0xff_b3b3b3
		};
		return colormap[index%colormap.length];
	}

	@Override
	public void paint(Graphics g) {
		super.paint(g);
		if(getWidth() < 1 || getHeight() < 1)
			return;
		g.setColor(getBackground());
		g.fillRect(0, 0, getWidth(), getHeight());
		paintDots((Graphics2D)g);
	}
	
	protected void paintDots(Graphics2D g2d) {
		// pointing upward y axis
		coordSysTransformGraphics(g2d);
		// scaling and translation so that content will be visible
		Rectangle2D view = getOptimalView();
		double translateX = view.getX();
		double translateY = view.getY();
		double scaleX = getWidth()/view.getWidth();
		double scaleY = getHeight()/view.getHeight();
		// drawing
		int k=0;
		for(double[][] points : pointSets) {
			int c = getColor(k);
			c = (c & 0x00_ffffff) | (pointOpacity << 24);
			g2d.setColor(new Color(c, true));
			for(int i=0; i<points.length; i++) {
				double x=getAt(points[i],0);
				double y=getAt(points[i],1);
				x-=translateX;
				y-=translateY;
				x*=scaleX;
				y*=scaleY;
				g2d.fillRect((int)x-1, (int)y-1, 2, 2);
			}
			k++;
		}
	}
	
	protected Rectangle2D getOptimalView() {
		Rectangle2D contentBounds = getContentBounds();
		double contentAspect = contentBounds.getWidth()/contentBounds.getHeight();
		double viewpAspect = getWidth()*1.0/getHeight();
		double w,h,x,y;
		if(viewpAspect < contentAspect) { // taller viewport
			w = contentBounds.getWidth();
			h = w/viewpAspect;
		} else { // wider viewport
			h = contentBounds.getHeight();
			w = h*viewpAspect;
		}
		x = contentBounds.getMinX()-(w-contentBounds.getWidth())/2;
		y = contentBounds.getMinY()-(h-contentBounds.getHeight())/2;
		return new Rectangle2D.Double(x, y, w, h);
	}
	
	void coordSysTransformGraphics(Graphics2D g2d) {
		g2d.translate(0, getHeight());
		g2d.scale(1.0, -1.0);
	}
	
	static double getAt(double[] arr, int idx) {
		return arr[idx%arr.length];
	}
	
	protected Rectangle2D calcContentBounds() {
		if(pointSets.isEmpty())
			return new Rectangle2D.Double(0, 0, 1, 1);
		double minX = pointSets.stream().flatMap(points->Arrays.stream(points))
				.mapToDouble(p->getAt(p,0)).min().getAsDouble();
		double maxX = pointSets.stream().flatMap(points->Arrays.stream(points))
				.mapToDouble(p->getAt(p,0)).max().getAsDouble();
		double minY = pointSets.stream().flatMap(points->Arrays.stream(points))
				.mapToDouble(p->getAt(p,1)).min().getAsDouble();
		double maxY = pointSets.stream().flatMap(points->Arrays.stream(points))
				.mapToDouble(p->getAt(p,1)).max().getAsDouble();
		return new Rectangle2D.Double(minX,minY, maxX-minX, maxY-minY);
	}
	
	public JFrame display(String title) {
		JFrame frame = new JFrame();
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.getContentPane().setLayout(new BorderLayout());
		frame.setTitle(title);
		frame.setMinimumSize(new Dimension(400, 400));
		
		frame.getContentPane().add(this, BorderLayout.CENTER);
		
		SwingUtilities.invokeLater(()->{
			frame.pack();
			frame.setVisible(true);
		});
		
		return frame;
	}
}
