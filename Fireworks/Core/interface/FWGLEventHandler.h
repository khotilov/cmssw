#ifndef Fireworks_Core_FWGLEventHandler_h
#define Fireworks_Core_FWGLEventHandler_h

#include "TGLEventHandler.h"
#include <sigc++/signal.h>

class TGWindow;

class FWGLEventHandler : public TGLEventHandler {
public:
   FWGLEventHandler(const char *name, TGWindow *w, TObject *obj, const char *title="");
   virtual ~FWGLEventHandler() {}

   virtual Bool_t HandleButton(Event_t * event);

   sigc::signal<void,Int_t,Int_t> openSelectedModelContextMenu_;

private:
   FWGLEventHandler(const FWGLEventHandler&); // stop default
   const FWGLEventHandler& operator=(const FWGLEventHandler&); // stop default
};

#endif
