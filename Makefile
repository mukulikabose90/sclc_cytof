IMAGE_NAME = sclc_cytof
WORKDIR = $(shell pwd)

build:
	docker build -t $(IMAGE_NAME) .

reproduce:
	docker run --rm \
		-v $(WORKDIR)/source:/project/source \
		-v $(WORKDIR)/data:/project/data \
		-v $(WORKDIR)/figures:/project/figures \
		$(IMAGE_NAME)

shell:
	docker run -it --rm \
		-v $(WORKDIR)/source:/project/source \
		-v $(WORKDIR)/data:/project/data \
		-v $(WORKDIR)/figures:/project/figures \
		$(IMAGE_NAME) bash

clean:
	docker rmi -f $(IMAGE_NAME)
