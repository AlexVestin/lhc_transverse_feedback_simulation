/* Note: this Google copyright notice only applies to the original file, which has large sections copy-pasted here. my changes are under CC0 (public domain).

 * Copyright 2015 Google Inc.
 *
 * Use of this source code is governed by a BSD-style license that can be
 * found in the LICENSE file.
 */

/*
The official instructions don't work well. These alternative instructions are intended to be the shortest path to get a minimal setup running.
The Linux steps were run through successfully on October 2019.
The Windows steps are known to be broken; the broken part is Step 7. The Include and Library directories should be tweaked.

This was made by copy-pasting and fixing two sources: https://github.com/google/skia/tree/master/experimental/GLFWTest and https://gist.github.com/zester/5163313
Don't bother trying these two sources; neither of them works.


step 1: install glfw (on Linux, "sudo apt install libglfw3-dev" will get you an acceptable (and outdated) version. on Visual Studio 2017, you must build glfw from source, contrary to Internet claims that glfw's VS2015 pre-compiled version works.)
step 2: follow the Setting Up section at http://commondatastorage.googleapis.com/chrome-infra-docs/flat/depot_tools/docs/html/depot_tools_tutorial.html#_setting_up
step 3: if you're in Windows, you will need a copy of bash; cmd.exe will fail in a later step. on my system, a copy of bash came with my installation of Git for Windows.
step 4: follow https://skia.org/user/download, using the "Clone the Skia repository" section only. Use Bash, even if you're on Windows. the Windows check "where python" is useful because sometimes python ends up in stupid places for stupid reasons
step 5: go to https://skia.org/user/build and look at the instructions, but don't follow them.
Move forward to either Windows step 6 or Linux step 6.


Windows step 6, Visual Studio 2017:
here is where bash is required, because cmd.exe doesn't allow single quotes, which are necessary to give the VC path. the various skia_use_foo commands are necessary to stop VS from erroring out when the headers are missing
run these two commands, replacing "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC" with your own VC directory:
gn gen out/Static --args='is_official_build=true win_vc="C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC" skia_use_libpng=false skia_use_zlib=false skia_use_libwebp=false skia_enable_pdf=false skia_use_libjpeg_turbo=false skia_use_expat=false'
ninja -C out/Static

Windows step 7. Warning, this is outdated. The massive "FOLDER1\x" was changed to not be necessary, but I haven't tested the correct steps now:
add this file to a new VS project
append "FOLDER1\skia\include\core;FOLDER1\skia\include\gpu;FOLDER1\skia\include\config;FOLDER1\skia\include\utils;FOLDER2\glfw\include'" to the VC include directories of your project, where FOLDERX represents the directories you put them in.
	you must include all 4 skia folders because the files inside skia folders assume they see the other folders.
	if you're unfamiliar with how the include directory works, it's in Project->Properties, VC++ Directories, Include Directories.
append "FOLDER1\skia\out\Static;FOLDER2\glfw\src\Debug;" to Library Directories, again replacing FOLDERX with the true location. add "opengl32.lib;skia.lib;glfw3.lib;" to Linker->Input->Additional Dependencies
Set build mode to x64.
Build! This will produce a debug mode binary.
If in the future you want a release mode binary, you will need to re-build glfw in release mode, and change the glfw library folder to FOLDER2\glfw\src\Release;

Linux step 6, Ubuntu 19.04. October 18, 2019:
Run:
sudo apt install clang libjpeg-dev libicu-dev libwebp-dev
bin/gn gen out/Static --args='is_official_build=true cc="clang" cxx="clang++"'
ninja -C out/Static

Linux step 7:
download this file as "glfw_ship.cpp", and place it in the parent folder of the "skia" directory. (this just makes "-Iskia" in the right place)
g++ -g -std=c++1z glfw_ship.cpp -lskia -ldl -lpthread -ljpeg -lfreetype -lz -lpng -lglfw -lfontconfig -lwebp -lwebpmux -lwebpdemux -lGL -Iskia -Lskia/out/Static/
./a.out


eventually, you will want color-correct spaces, and there are 5 places below (Ctrl+F "enable correct color spaces"), where you should replace/uncomment lines to enable this.
warning: color-correct spaces don't work in VMWare, because mesa doesn't support it.
*/
#define SK_GL 1
#include <iostream>
#include <sstream>
#include <GLES3/gl3.h>
#include "GLFW/glfw3.h"
#include "include/gpu/GrBackendSurface.h"
#include "include/gpu/GrDirectContext.h"
#include "include/core/SkTextBlob.h"

#include "include/gpu/gl/GrGLInterface.h"
#include "include/core/SkCanvas.h"
#include "include/core/SkColorSpace.h"
#include "include/core/SkSurface.h"
#include <stdio.h>
#include <stdlib.h>

//uncomment the two lines below to enable correct color spaces
//#define GL_FRAMEBUFFER_SRGB 0x8DB9
//#define GL_SRGB8_ALPHA8 0x8C43
#define GR_GL_RGBA8                          0x8058
#define GR_GL_BGRA8                          0x93A1

sk_sp<GrDirectContext> dContext;
sk_sp<SkSurface> dSurface;
SkCanvas* dCanvas;
GLFWwindow* window;


double firstFrequency = -1;
const bool usefirstFrequency = true;
int counter = 0;
double avgError = 0, currentAvgError = 0;



int frequencyToBin(double frequency, int nfft, int sampleRate) {
    double fRes = sampleRate / (double)nfft;
    return frequency / fRes; 
}

double binToFrequency(int bin, int nfft, int sampleRate) {
    double fRes = sampleRate / (double)nfft;
    return fRes * bin; 
}

void error_callback(int error, const char* description) {
	fputs(description, stderr);
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
		glfwSetWindowShouldClose(window, GL_TRUE);
}

void init_skia(int w, int h) {

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glClearColor(0, 0, 0, 0);
    glClearStencil(0);
    glClear(GL_COLOR_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

    glEnable(GL_STENCIL_TEST);    


	GrContextOptions options;
	//options.fRequireDecodeDisableForSRGB = false; //was removed?
    auto interface = GrGLMakeNativeInterface();
    dContext = GrDirectContext::MakeGL(interface);

    dContext->resetContext(kRenderTarget_GrGLBackendState | kMisc_GrGLBackendState);


    GrGLFramebufferInfo info;
    info.fFBOID = 0;
    info.fFormat = GR_GL_RGBA8;

    GrGLint sampleCnt = 0;
    // glGetIntegerv(GL_SAMPLES, &sampleCnt);

    GrGLint stencil = 8;
    // TODO not valid
    // glGetIntegerv(GL_STENCIL_BITS, &stencil);

    std::cout << "Stencil: " << stencil << " sampleCnt: " << sampleCnt << std::endl;
    if(stencil == 0) {
        std::cerr << "Couldn't create stencil buffer" << std::endl;
        exit(1);
    }
	SkColorType colorType = kBGRA_8888_SkColorType;
	GrBackendRenderTarget target(w, h, sampleCnt, stencil, info);

    //SkSurfaceProps props;


	//(replace line below with this one to enable correct color spaces) sSurface = SkSurface::MakeFromBackendRenderTarget(sContext, backendRenderTarget, kBottomLeft_GrSurfaceOrigin, colorType, SkColorSpace::MakeSRGB(), nullptr).release();
	dSurface = SkSurface::MakeFromBackendRenderTarget(
        dContext.get(), 
        target,
        kBottomLeft_GrSurfaceOrigin, 
        colorType, 
        nullptr, 
        nullptr
    );

	if (dSurface == nullptr) {
        std::cout << "Surface not created" << std::endl;
        exit(1);
    }

}

void cleanup_skia() {
    // TODO cleanup
	//dSurface.delete();
	//delete dContext;
}

const int kWidth = 1200;
const int kHeight = 640;


void DrawString(const char* text, SkFont& font, SkPaint& paint, int x, int y) {
    auto blob = SkTextBlob::MakeFromText(text, strlen(text), font, SkTextEncoding::kUTF8);
    dCanvas->drawTextBlob(blob, x, y, paint);
}

void ClearCanvas() {
    dCanvas->clear( { 1.0, 1.0, 1.0, 1.0 });
    glfwPollEvents();
}

void DrawPoints(std::vector<double> data, double realFrequency, int nfft, int sampleRate, SkColor4f col) {
    SkPaint paint;

    dCanvas->save();
    // Move to middle of canvas
    dCanvas->translate(kWidth/2, kHeight/2);
    // freq res and stuff
    double frequencyResolution = sampleRate / (double)nfft;

    // Find max, index of max and mean
    int maxIndex = -1;
    double maxValue = -1000;
    double avg = 0;
    for(int i = 0; i < data.size(); i++) {
        avg += data[i];
        if(data[i] > maxValue) {
            maxValue = data[i];
            maxIndex = i;
        }
    }
    
    // Sizing of graph
    int size = data.size();
    

    if (firstFrequency == -1) {
        firstFrequency = maxIndex;
    }

    int frange = 20;

    double lf = realFrequency - frange;
    double uf = realFrequency + frange;

    int lower = std::max(frequencyToBin(lf, nfft, sampleRate), 0);
    int upper = std::min(frequencyToBin(uf, nfft, sampleRate), size - 1);
    int nPoints = upper - lower;

    double width = kWidth / 2.;


    // If scale for maxvalue
    // const double height = ((float)kHeight / 4.) / maxValue;
    double height = ((float)kHeight / 4.) * 0.25;


    // Box sizes
    int s = 4;
    int ms = 6;
    dCanvas->translate(-kWidth/4, 0);
    for(int i = 0; i < nPoints; i++) {
        int idx = i + lower; 
        double f = binToFrequency(idx, nfft, sampleRate);
        double pos = (f - lf) / (uf - lf);  

        double h = -data[idx] * height;


        if (idx == maxIndex) {
            paint.setColor4f({ 0, 1, 0, 1.0 });
            float x1 = pos*width - ms/2.;
            float y1 = h + ms/2.;
            float x2 = pos*width + ms/2.;
            float y2 = h - ms/2.;

            dCanvas->drawRect( { x1, y1, x2, y2 }, paint);
        }

        paint.setColor4f(col);

        float x1 = pos*width - s/2;
        float y1 = h + s/2;
        float x2 = pos*width + s/2;
        float y2 = h - s/2;
        dCanvas->drawRect( { x1, y1, x2, y2 }, paint);
    }   

    dCanvas->translate(kWidth/4, 0);

    paint.setColor4f({1.0, 0.0, 0.0, 0.5 }); 
    paint.setStyle(SkPaint::kStroke_Style);
    SkPath path;
    path.moveTo(0, 0);
    path.lineTo(0, -150);
    paint.setStrokeWidth(1);
    dCanvas->drawPath(path, paint);

    /*
    // Draw text stuff
    SkPaint textPaint;
    SkFont font;
    const int textPosX = - ((float)kWidth / 4.);
    const int textPosY = 20;

    const int fontSize = 16;
    const int paddedFontSize = 18;
    font.setSize(fontSize);

    std::string infoStr = "Index: " + 
        std::to_string(maxIndex) + 
        " freq resolution: " + 
        std::to_string(frequencyResolution) +
        " max value: " +
        std::to_string(maxValue);
    DrawString(infoStr.c_str(), font, textPaint, textPosX, textPosY);

    std::string maxFreqStr = "Found max frequency: " + std::to_string(frequencyResolution * maxIndex);
    DrawString(maxFreqStr.c_str(), font, textPaint, textPosX, textPosY + paddedFontSize * 1);

    
    std::string realFreqStr = "Real frequency: " + std::to_string(realFrequency);
    DrawString(realFreqStr.c_str(), font, textPaint, textPosX, textPosY + paddedFontSize * 2);


    textPaint.setColor4f({ 1.0, 0.0, 0.0, 1.0 });
    counter++;

    double error = std::abs(frequencyResolution * maxIndex - realFrequency);
    avgError += error;

    std::string errorStr = "Error: " + std::to_string(error) + " Average error: " + std::to_string(avgError / (double)counter);    
    DrawString(errorStr.c_str(), font, textPaint, textPosX, textPosY + paddedFontSize * 3);
    */
    dCanvas->restore();

}

bool FlushCanvas() {
    dContext->flush();
    glfwSwapBuffers(window);
    return !glfwWindowShouldClose(window); 
}


void CloseWindow() {
    cleanup_skia();
	glfwDestroyWindow(window);
	glfwTerminate();
	exit(EXIT_SUCCESS);
}

int InitWindow(void) {
	glfwSetErrorCallback(error_callback);
	if (!glfwInit()) {
		exit(EXIT_FAILURE);
	}

	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 5);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	//(uncomment to enable correct color spaces) 
    glfwWindowHint( GLFW_STENCIL_BITS, 8 );
    glfwWindowHint( GLFW_DEPTH_BITS , 24 );

	window = glfwCreateWindow(kWidth, kHeight, "Simple example", NULL, NULL);
	if (!window) {
		glfwTerminate();
		exit(EXIT_FAILURE);
	}
	glfwMakeContextCurrent(window);
	//(uncomment to enable correct color spaces) glEnable(GL_FRAMEBUFFER_SRGB);
	init_skia(kWidth, kHeight);


	glfwSwapInterval(1);
	glfwSetKeyCallback(window, key_callback);

	// Draw to the surface via its SkCanvas.
	dCanvas = dSurface->getCanvas(); // We don't manage this pointer's lifetime.
    return 0;
}