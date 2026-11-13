/*        L O D _ P R O G R E S S _ O V E R L A Y . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "BObol/BLodProgressOverlay.h"
#include "BObol/BViewController.h"

#include <Inventor/annex/HUD/nodekits/SoHUDKit.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/elements/SoCacheElement.h>
#include <Inventor/elements/SoFontSizeElement.h>
#include <Inventor/elements/SoLazyElement.h>
#include <Inventor/elements/SoModelMatrixElement.h>
#include <Inventor/elements/SoViewportRegionElement.h>
#include <Inventor/fields/SoSFBool.h>
#include <Inventor/fields/SoSFColor.h>
#include <Inventor/fields/SoSFString.h>
#include <Inventor/misc/SoState.h>
#include <Inventor/nodes/SoNode.h>
#include <Inventor/nodes/SoText2.h>
#include <Inventor/SbViewportRegion.h>
#include <Inventor/system/SoGLDispatch.h>
#include <Inventor/system/gl.h>
#include <Inventor/tools/SbModernUtils.h>

#include <algorithm>
#include <chrono>

namespace {

constexpr float CARD_WIDTH = 420.0f;
constexpr float CARD_HEIGHT = 30.0f;
constexpr float CARD_MARGIN = 8.0f;
constexpr float CARD_PADDING = 10.0f;
constexpr float TEXT_FONT_SIZE = 11.0f;
constexpr float SPINNER_RADIUS = 6.0f;
constexpr float SPINNER_DOT_SIZE = 2.5f;
constexpr int SPINNER_DOT_COUNT = 8;
constexpr int SPINNER_STEP_MILLISECONDS = 90;

class SoBRLLodProgressCard : public SoNode {
    typedef SoNode inherited;

    SO_NODE_HEADER(SoBRLLodProgressCard);

public:
    SoSFString title;
    SoSFString detail;
    SoSFColor color;
    SoSFBool terminal;

    SoBRLLodProgressCard(void) : textNode(NULL)
    {
	SO_NODE_CONSTRUCTOR(SoBRLLodProgressCard);
	SO_NODE_ADD_FIELD(title, (""));
	SO_NODE_ADD_FIELD(detail, (""));
	SO_NODE_ADD_FIELD(color, (1.0f, 0.75f, 0.28f));
	SO_NODE_ADD_FIELD(terminal, (FALSE));
    }

    static void initClass(void)
    {
	SO_NODE_INIT_CLASS(SoBRLLodProgressCard, SoNode, "Node");
    }

    void GLRender(SoGLRenderAction *action) override
    {
	SoState *state = action->getState();
	if (!this->terminal.getValue())
	    SoCacheElement::invalidate(state);
	const SbVec2s viewport =
	    SoViewportRegionElement::get(state).getViewportSizePixels();
	const float availableWidth =
	    std::max(0.0f, static_cast<float>(viewport[0]) - 2.0f * CARD_MARGIN);
	const float width = std::min(CARD_WIDTH, availableWidth);
	if (width <= 2.0f * CARD_PADDING)
	    return;

	const float left = CARD_MARGIN;
	const float bottom = static_cast<float>(viewport[1]) - CARD_MARGIN -
	    CARD_HEIGHT;
	const float right = left + width;
	const float top = bottom + CARD_HEIGHT;
	const SbColor accent = this->color.getValue();
	const SoGLContext *gl = sogl_glue_from_state(state);

	SoGLContext_glPushAttrib(gl, GL_ENABLE_BIT | GL_COLOR_BUFFER_BIT |
	    GL_CURRENT_BIT | GL_LINE_BIT);
	SoGLContext_glDisable(gl, GL_LIGHTING);
	SoGLContext_glDisable(gl, GL_TEXTURE_2D);
	SoGLContext_glEnable(gl, GL_BLEND);
	SoGLContext_glBlendFunc(gl, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	SoGLContext_glColor4f(gl, 0.055f, 0.071f, 0.090f, 0.86f);
	drawQuad(gl, left, bottom, right, top);
	SoGLContext_glColor4f(gl, accent[0], accent[1], accent[2], 0.90f);
	drawQuad(gl, left, bottom, left + 2.0f, top);

	const float centerX = left + CARD_PADDING + SPINNER_RADIUS;
	const float centerY = bottom + CARD_HEIGHT * 0.5f;
	if (this->terminal.getValue()) {
	    drawQuad(gl, centerX - 4.0f, centerY - 4.0f,
		centerX + 4.0f, centerY + 4.0f);
	} else {
	    static const float offsets[SPINNER_DOT_COUNT][2] = {
		{0.0f, 1.0f}, {0.707f, 0.707f}, {1.0f, 0.0f},
		{0.707f, -0.707f}, {0.0f, -1.0f}, {-0.707f, -0.707f},
		{-1.0f, 0.0f}, {-0.707f, 0.707f}
	    };
	    const auto elapsed = std::chrono::duration_cast<
		std::chrono::milliseconds>(
		    std::chrono::steady_clock::now().time_since_epoch()).count();
	    const int activeDot = static_cast<int>(
		(elapsed / SPINNER_STEP_MILLISECONDS) % SPINNER_DOT_COUNT);
	    for (int i = 0; i < SPINNER_DOT_COUNT; ++i) {
		const int age = (activeDot - i + SPINNER_DOT_COUNT) %
		    SPINNER_DOT_COUNT;
		const float alpha = std::max(0.18f, 1.0f - 0.12f * age);
		const float dotX = centerX + offsets[i][0] * SPINNER_RADIUS;
		const float dotY = centerY + offsets[i][1] * SPINNER_RADIUS;
		SoGLContext_glColor4f(gl, accent[0], accent[1], accent[2],
		    alpha);
		drawQuad(gl, dotX - SPINNER_DOT_SIZE * 0.5f,
		    dotY - SPINNER_DOT_SIZE * 0.5f,
		    dotX + SPINNER_DOT_SIZE * 0.5f,
		    dotY + SPINNER_DOT_SIZE * 0.5f);
	    }
	}
	SoGLContext_glPopAttrib(gl);

	SbString line = this->title.getValue();
	if (this->detail.getValue().getLength() > 0) {
	    line += "  -  ";
	    line += this->detail.getValue();
	}
	renderText(action, this->textNode, line,
	    SbVec2f(centerX + SPINNER_RADIUS + CARD_PADDING,
		bottom + 10.0f), TEXT_FONT_SIZE,
	    SbColor(0.90f, 0.92f, 0.95f));
    }

protected:
    ~SoBRLLodProgressCard(void) override
    {
	if (this->textNode)
	    this->textNode->unref();
    }

private:
    static void drawQuad(const SoGLContext *gl, float left, float bottom,
	float right, float top)
    {
	SoGLContext_glBegin(gl, GL_QUADS);
	SoGLContext_glVertex2f(gl, left, bottom);
	SoGLContext_glVertex2f(gl, right, bottom);
	SoGLContext_glVertex2f(gl, right, top);
	SoGLContext_glVertex2f(gl, left, top);
	SoGLContext_glEnd(gl);
    }

    void renderText(SoGLRenderAction *action, SoText2 *&textNode,
	const SbString &text, const SbVec2f &position, float fontSize,
	const SbColor &textColor)
    {
	SoState *state = action->getState();
	state->push();
	SoModelMatrixElement::makeIdentity(state, this);
	SoModelMatrixElement::translateBy(state, this,
	    SbVec3f(position[0], position[1], 0.0f));
	SoFontSizeElement::set(state, this, fontSize);
	SoColorPacker colorPacker;
	SoLazyElement::setDiffuse(state, this, 1, &textColor, &colorPacker);
	if (!textNode) {
	    textNode = new SoText2;
	    textNode->ref();
	}
	textNode->string.setValue(text);
	textNode->GLRender(action);
	state->pop();
    }

    SoText2 *textNode;
};

SO_NODE_SOURCE(SoBRLLodProgressCard);

}

SO_NODE_SOURCE(SoBRLLodProgressOverlay);

SoBRLLodProgressOverlay::SoBRLLodProgressOverlay(void)
{
    SO_NODE_CONSTRUCTOR(SoBRLLodProgressOverlay);

    SO_NODE_ADD_FIELD(title, (""));
    SO_NODE_ADD_FIELD(detail, (""));
    SO_NODE_ADD_FIELD(color, (1.0f, 0.75f, 0.28f));
    SO_NODE_ADD_FIELD(visible, (TRUE));
    SO_NODE_ADD_FIELD(terminal, (FALSE));
    SO_NODE_ADD_FIELD(terminalReady, (FALSE));
}

SoBRLLodProgressOverlay::~SoBRLLodProgressOverlay(void)
{
}

void
SoBRLLodProgressOverlay::initClass(void)
{
    SoBRLLodProgressCard::initClass();
    SO_NODE_INIT_CLASS(SoBRLLodProgressOverlay, SoSeparator, "Separator");
}

void
SoBRLLodProgressOverlay::setStatus(
    const BObolLodProgressPresentationStatus &status)
{
    this->title = status.title;
    this->detail = status.detail;
    this->color = status.color;
    this->visible = status.visible;
    this->terminal = status.terminal;
    this->terminalReady = status.terminalReady;
}

SoHUDKit *
SoBRLLodProgressOverlay::rebuildGeometry(void)
{
    if (!this->visible.getValue() || this->title.getValue().getLength() == 0) {
	auto replacement = this->getChildren()->prepareReplacement({});
	replacement->commit();
	replacement->notify();
	return NULL;
    }

    SbModernUtils::SoNodeRef hudOwner(new SoHUDKit);
    auto *hud = static_cast<SoHUDKit *>(hudOwner.get());
    SbModernUtils::SoNodeRef cardOwner(new SoBRLLodProgressCard);
    auto *card = static_cast<SoBRLLodProgressCard *>(cardOwner.get());
    card->title = this->title.getValue();
    card->detail = this->detail.getValue();
    card->color = this->color.getValue();
    card->terminal = this->terminal.getValue();
    hud->addWidget(card);

    auto replacement = this->getChildren()->prepareReplacement({hud});
    replacement->commit();
    replacement->notify();
    return hud;
}

SoHUDKit *
SoBRLLodProgressOverlay::getHUDKit(void) const
{
    for (int i = 0; i < this->getNumChildren(); i++) {
	SoNode *node = this->getChild(i);
	if (node && node->isOfType(SoHUDKit::getClassTypeId()))
	    return static_cast<SoHUDKit *>(node);
    }
    return NULL;
}

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
