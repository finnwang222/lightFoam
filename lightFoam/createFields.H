

            volScalarField  IDefault
            (
                IOobject
                (
                    "IDefault",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
               mesh
            );




    autoPtr<light::lightModel> lightM
    (
       light::lightModel::New(IDefault)
    );

