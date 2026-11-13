/*                        I S S T G L . H
 * BRL-CAD
 *
 * Copyright (c) 2021-2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this file; see the file named COPYING for more
 * information.
 */
/** @file isstgl.h
 *
 * QWidget image presenter encoding the information specific to the
 * TIE raytracing view.
 *
 * TODO:  Look at f_knob, knob_rot, mged_vrot_syz, mged_rot, etc.
 * to determine how MGED is doing its mouse x,y delta to view
 * manipulation translation logic.
 *
 */

#ifndef ISSTGL_H
#define ISSTGL_H

#include <QImage>
#include <QWidget>
#include <QThread>

#include <atomic>

extern "C" {
#include "rt/tie.h"
#include "librender/camera.h"
}

class isstView;

class TIERenderer : public QObject
{
    Q_OBJECT
    public:
	TIERenderer();
	~TIERenderer();

	void prepareExit() { m_exiting.store(true); }

    signals:
	void imageReady(const QImage &image);

    public slots:
	void render();
	void res_incr();
	void res_decr();
	void setResolution(int resolution);
	void setSize(int width, int height);
	void setTie(struct tie_s *in_tie);
	void clearTie();

    public:
	struct tie_s *tie = NULL; // From parent app
	struct render_camera_s camera;
	vect_t camera_pos_init;
	vect_t camera_focus_init;
	bool changed = true;

    private:
	bool resize();

	struct camera_tile_s tile;
	tienet_buffer_t buffer_image;
	int resolution = 20;
	int resolution_factor = 0;
	int viewport_width = 0;
	int viewport_height = 0;
	std::atomic_bool m_exiting{false};
};

class isstView : public QWidget
{
    Q_OBJECT

    public:
	explicit isstView(QWidget *parent = nullptr);
	~isstView();

	void set_tie(struct tie_s *in_tie);
	void clear_tie();
	void setPreviewResolution(int resolution);

	void save_image();

    protected:
	void paintEvent(QPaintEvent *event) override;
	void resizeEvent(QResizeEvent *event) override;

	void keyPressEvent(QKeyEvent *k) override;

    signals:
      void resolutionIncreased();
      void resolutionDecreased();
      void resolutionRequested(int resolution);
      void sceneChanged(struct tie_s *in_tie);
      void sizeChanged(int width, int height);
      void imagePresented();

    public slots:
      void setImage(const QImage &image);

    private:
	QThread *m_thread;
	TIERenderer *m_renderer;
	QImage m_image;
};

#endif /* ISSTGL_H */

// Local Variables:
// tab-width: 8
// mode: C++
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
