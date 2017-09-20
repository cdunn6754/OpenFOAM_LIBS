09-16-17
All of the DAEM's do have a purpose. We can not just make a new devol model type becase the original base class DevolatilisationModel, does not have support for integrated rates. So I made a whole new base class and the corresponding additional models (Single, Constant) in addition to the one we really wanted DAEM. I added the DAEM to all models and file names right before the final Devolatilisation to distiguish the two groups (Daem vs original) of classes. 

09-20-17
Update, i had to change the typenames to the originals but left the actual clas names prefixed with Daem. That way the new daem parcels/clouds can call the same name in coalCloudProperties but get the new classes. 

Not entirely sure if changing all of the models (the class names of the models) was necessary.
