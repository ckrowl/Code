#pragma once

#include "matrix4.h"

struct RigidBody
{
	RigTForm rtf;
	Matrix4 scale;
	RigidBody **children;
	int numOfChildren;
	Cvec3 color;
	Geometry *geom;
	bool isVisible;
	string name;

	RigidBody()
	{
		rtf = RigTForm();
		scale = Matrix4();
		children = NULL;
		numOfChildren = 0;
		color = Cvec3(.5,.5,.5);
		geom = NULL;
		isVisible = true;
	}

	~RigidBody()
	{
		if(children)
		{
			for(int i = 0; i < numOfChildren; ++i)
				delete children[i];
			delete [] children;
		}
	}

	RigidBody(RigTForm rtf_, Matrix4 scale_, RigidBody **children_, Geometry *geom_, Cvec3 color_)
	{
		/* PURPOSE:		 
			RECEIVES:	 
							
			RETURNS:		
		*/

		rtf = rtf_;
		scale = scale_;
		children = children_;
		numOfChildren = 0;
		geom = geom_;
		color = color_;
		isVisible = true;
	}

	void drawRigidBody(const ShaderState& curSS, RigTForm invEyeRbt)
	{
		RigTForm respectFrame = invEyeRbt;// * rtf;
		draw(curSS, respectFrame, Matrix4());
		
		//2nd Method Rotations are scaled correctly Other scaling problems
		rtf = rtf * RigTForm(Cvec3(3,0,0));
		Matrix4 respectFrame2 = RigTForm::makeTRmatrix(invEyeRbt,scale);
		draw(curSS, respectFrame2);
		rtf = rtf * RigTForm(Cvec3(-3,0,0));
	}

	void draw(const ShaderState& curSS, Matrix4 respectFrame_)
	{
		safe_glUniform3f(curSS.h_uColor, color[0], color[1], color[2]);
			
		//Draw parent
		Matrix4 respectFrame = RigTForm::makeTRmatrix(rtf,scale);
		respectFrame *= respectFrame_;
		Matrix4 MVM = respectFrame;

		if (isVisible)
		{
			if (geom != NULL)
				geom->draw(curSS, MVM);
		}

		//Draw Children
		for (int i = 0; i < numOfChildren; i++)
		{
			children[i]->draw(curSS, respectFrame);
		}
		
	}

	void draw(const ShaderState& curSS, RigTForm respectFrame_, Matrix4 respectScale_)
	{
		safe_glUniform3f(curSS.h_uColor, color[0], color[1], color[2]);
			
		//Draw parent
		this;
	
		//scale correct but not translated correctly; moving one object scale moves the rest;
		RigTForm respectFrame = respectFrame_ * rtf;
		Matrix4 respectScale = respectScale_ * scale;
		Matrix4 MVM = RigTForm::makeTRmatrix(respectFrame,respectScale);
		
		/*
		//Positioning doesn't change after scales; Moving one object doesn't translate children during setup
		RigTForm respectFrame = respectFrame_ * rtf;
		Matrix4 respectScale = respectScale_ * scale;
		Matrix4 temp1 = RigTForm::makeTRmatrix(respectFrame_) * respectScale_;
		Matrix4 temp2 = RigTForm::makeTRmatrix(rtf) * scale;
		Matrix4 MVM = temp1 * temp2;
		*/

		if (isVisible)
		{
			if (geom != NULL)
				geom->draw(curSS, MVM);
		}

		//Draw Children
		for (int i = 0; i < numOfChildren; i++)
		{
			children[i]->draw(curSS, respectFrame, respectScale);
		}
		
	}
};