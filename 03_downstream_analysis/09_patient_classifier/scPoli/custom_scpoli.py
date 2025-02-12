from scarches.trainers.scpoli.trainer import scPoliTrainer
from scarches.models.scpoli.scpoli_model import scPoli
import numpy as np

from scarches.dataset.scpoli.anndata import MultiConditionAnnotatedDataset

class CustomScPoliTrainer(scPoliTrainer):

    def __init__(self, *args, w_logger, prefix, **kwargs):
        super().__init__(*args, **kwargs)
        self.w_logger = w_logger
        self.steps = -1
        self.prefix = prefix
        self.w_logger.define_metric("global_step")
        self.w_logger.define_metric("*", step_metric="global_step", step_sync=True)

    def on_iteration(self, batch_data):
        super().on_iteration(batch_data)
        self.steps += 1
        key = self.prefix + "train_loss"
        self.w_logger.log({key: self.current_loss, 'global_step':self.steps})

    def on_epoch_end(self):
        super().on_epoch_end()
        
        self.w_logger.log({self.prefix+'epoch': self.epoch, 'global_step':self.steps})

        for metric in self.iter_logs:
            if metric.startswith('val_'):
                self.w_logger.log({self.prefix + metric: self.iter_logs[metric][-1], 'global_step':self.steps})
            else:
                self.w_logger.log({self.prefix + "epoch_" + metric: self.iter_logs[metric][-1], 'global_step':self.steps})
                
    def loss(self, total_batch=None):
        latent, recon_loss, kl_loss, mmd_loss = self.model(**total_batch)
        self.iter_logs["recon_loss"].append(recon_loss.item())
        self.iter_logs["kl_loss"].append(kl_loss.item())
        self.iter_logs["mmd_loss"].append(mmd_loss.item())

        #calculate classifier loss for labeled/unlabeled data
        label_categories = total_batch["labeled"].unique().tolist()
        unweighted_prototype_loss = torch.tensor(0.0, device=self.device)
        unlabeled_loss = torch.tensor(0.0, device=self.device)
        labeled_loss = torch.tensor(0.0, device=self.device)
        if self.epoch >= self.pretraining_epochs:
            #calculate prototype loss for all data
            if self.prototypes_unlabeled is not None:
                unlabeled_loss, _ = self.prototype_unlabeled_loss(
                    latent,
                    torch.stack(self.prototypes_unlabeled).squeeze(),
                )
                unweighted_prototype_loss = (
                    unweighted_prototype_loss + self.unlabeled_weight * unlabeled_loss
                )

            # Calculate prototype loss for labeled data
            if (self.any_labeled_data is True) and (self.prototype_training is True):
                labeled_loss = self.prototype_labeled_loss(
                    latent[torch.where(total_batch["labeled"] == 1)[0], :],
                    self.prototypes_labeled,
                    total_batch["celltypes"][
                        torch.where(total_batch["labeled"] == 1)[0], :
                    ],
                )
                unweighted_prototype_loss = unweighted_prototype_loss + labeled_loss

        # Loss addition and Logs
        prototype_loss = self.eta * unweighted_prototype_loss
        cvae_loss = recon_loss + self.calc_alpha_coeff() * kl_loss + mmd_loss
        loss = cvae_loss + prototype_loss
        self.iter_logs["loss"].append(loss.item())
        self.iter_logs["unweighted_loss"].append(
            recon_loss.item()
            + kl_loss.item()
            + mmd_loss.item()
            + unweighted_prototype_loss.item()
        )
        self.iter_logs["cvae_loss"].append(cvae_loss.item())
        if self.epoch >= self.pretraining_epochs:
            self.iter_logs["prototype_loss"].append(prototype_loss.item())
            if 0 in label_categories or self.model.unknown_ct_names is not None:
                self.iter_logs["unlabeled_loss"].append(unlabeled_loss.item())
            if 1 in label_categories:
                self.iter_logs["labeled_loss"].append(labeled_loss.item())
        return loss

class CustomScPoli(scPoli):
    def train(
        self,
        n_epochs: int = 100,
        pretraining_epochs=None,
        eta: float = 1,
        lr: float = 1e-3,
        eps: float = 0.01,
        alpha_epoch_anneal = 1e2,
        reload_best: bool = False,
        prototype_training = True,
        unlabeled_prototype_training = True,
        w_logger = None,
        prefix = None,
        **kwargs,
    ):
        """Train the model.

        Parameters
        ----------
        n_epochs
             Number of epochs for training the model.
        lr
             Learning rate for training the model.
        eps
             torch.optim.Adam eps parameter
        kwargs
             kwargs for the scPoli trainer.
        """
        self.prototype_training_ = prototype_training
        self.unlabeled_prototype_training_ = unlabeled_prototype_training
        if self.cell_type_keys_ is None:
            pretraining_epochs = n_epochs
            self.prototype_training_ = False
            print("The model is being trained without using prototypes.")
        elif pretraining_epochs is None:
            pretraining_epochs = int(np.floor(n_epochs * 0.9))


        self.trainer = CustomScPoliTrainer(
            self.model,
            self.adata,
            labeled_indices=self.labeled_indices_,
            pretraining_epochs=pretraining_epochs,
            condition_keys=self.condition_keys_,
            cell_type_keys=self.cell_type_keys_,
            reload_best=reload_best,
            prototype_training=self.prototype_training_,
            unlabeled_prototype_training=self.unlabeled_prototype_training_,
            eta=eta,
            alpha_epoch_anneal=alpha_epoch_anneal,
            w_logger=w_logger,
            prefix = prefix,
            **kwargs,
        )
        self.trainer.train(n_epochs, lr, eps)
        self.is_trained_ = True
        self.prototypes_labeled_ = self.model.prototypes_labeled
        self.prototypes_unlabeled_ = self.model.prototypes_unlabeled
