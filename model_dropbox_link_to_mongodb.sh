# Run using `bash model_dropbox_link_to_mongodb.sh <MODEL_DIRECTORY_BASENAME>`

~/Dropbox-Uploader/dropbox_uploader.sh upload models/uvvis/production/$1.tar.gz /
SHAREABLE_LINK=$(~/Dropbox-Uploader/dropbox_uploader.sh share $1.tar.gz)
python model_dropbox_link_to_mongodb.py $SHAREABLE_LINK uv_vis
